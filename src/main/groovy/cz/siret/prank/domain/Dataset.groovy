package cz.siret.prank.domain

import com.google.common.base.CharMatcher
import com.google.common.base.Splitter
import cz.siret.prank.domain.loaders.ConcavityLoader
import cz.siret.prank.domain.loaders.FPockeLoader
import cz.siret.prank.domain.loaders.PredictionLoader
import cz.siret.prank.domain.loaders.SiteHoundLoader
import cz.siret.prank.features.api.ProcessedItemContext
import cz.siret.prank.program.PrankException
import cz.siret.prank.program.ThreadPoolFactory
import cz.siret.prank.program.params.Parametrized
import cz.siret.prank.utils.Futils
import cz.siret.prank.utils.StrUtils
import groovy.util.logging.Slf4j
import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger

/**
 * Dataset represents a list of items (usually proteins) to be processed by the program.
 * Multi-column format with declared variable header allows to specify complementary data.
 *
 * see distro/test_data/readme.txt for dataset format specification
 */
@Slf4j
class Dataset implements Parametrized {

    static final Splitter SPLITTER = Splitter.on(CharMatcher.whitespace()).trimResults().omitEmptyStrings()

    /*
     * dataset parameter names
     * Dataset parameters an be defined in dataset file as PARAM.<PARAM_NAME>=<value>
     */
    static final String PARAM_PREDICTION_METHOD = "PREDICTION_METHOD"
    static final String PARAM_LIGANDS_SEPARATED_BY_TER = "LIGANDS_SEPARATED_BY_TER"

    /*
     * dataset column names
     */
    static final String COLUMN_PROTEIN = "protein"
    static final String COLUMN_PREDICTION = "prediction"
    static final String COLUMN_LIGAND_CODES = "ligand_codes"

    static final List<String> DEFAULT_HEADER = [ COLUMN_PROTEIN ]


    /**
     * Contains file names, (optionally) ligand codes and cached structures.
     */
    class Item {
        Map<String, String> columnValues
        String proteinFile  // liganated/unliganated protein for predictions
        String pocketPredictionFile // nullable

        String label
        PredictionPair cachedPair

        private Item(String proteinFile, String predictionFile, Map<String, String> columnValues) {
            this.columnValues = columnValues
            this.proteinFile = proteinFile
            this.pocketPredictionFile = predictionFile
            this.label = Futils.shortName( pocketPredictionFile ?: proteinFile )
        }

        Prediction getPrediction() {
            // when running 'prank rescore' on a dataset with one column prediction is in proteinFile
            String file = pocketPredictionFile ?: proteinFile
            return getLoader(this).loadPredictionWithoutProtein(file)
        }

        PredictionPair getPredictionPair() {
            if (cached) {
                if (cachedPair==null) {
                    cachedPair = loadPredictionPair()
                    log.info "caching structures in dataset item [$label]"
                }
                return cachedPair
            } else {
                return loadPredictionPair()
            }
        }

        // for one column datasets
        Protein getProtein() {
            getPredictionPair().liganatedProtein
        }

        PredictionPair loadPredictionPair() {
            getLoader(this).loadPredictionPair(proteinFile, pocketPredictionFile)
        }

        /**
         * explicitely specified ligand codes
         * @return null if column is not defined
         */
        Set<String> getLigandCodes() {
            if (!columnValues.containsKey(COLUMN_LIGAND_CODES)) {
                null
            } else {
                Splitter.on(",").split(columnValues[COLUMN_LIGAND_CODES]).toList()
            }
        }

        ProcessedItemContext getContext() {
            new ProcessedItemContext(columnValues)
        }
    }

    Item createNewItem(String proteinFile, String predictionFile, Map<String, String> columnValues) {
        return new Item(proteinFile, predictionFile, columnValues)
    }

    Item createNewItem(Map<String, String> columnValues) {
        String proteinFile = dir + "/" + columnValues.get(COLUMN_PROTEIN)
        String predictionFile = null
        if (header.contains(COLUMN_PREDICTION)) {
            predictionFile = dir + "/" + columnValues.get(COLUMN_PREDICTION)
        }

        return createNewItem(predictionFile, predictionFile, columnValues)
    }

    static final class Fold {
        int num
        Dataset trainset
        Dataset evalset
        Fold(int num, Dataset trainset, Dataset evalset) {
            this.num = num
            this.trainset = trainset
            this.evalset = evalset
        }
    }

    interface Processor {

        abstract void processItem(Item item)
    }

    /**
     * summary of a dataset processing run
     */
    static class Result {
        List<Item> errorItems = Collections.synchronizedList(new ArrayList<Item>())

        boolean hasErrors() {
            errorItems.size() > 0
        }

        int getErrorCount() {
            errorItems.size()
        }
    }

//===========================================================================================================//

    String name
    String dir
    Map<String, String> attributes = new HashMap<>()
    List<String> header = DEFAULT_HEADER
    List<Item> items = new ArrayList<>()
    boolean cached = false

//===========================================================================================================//

    Dataset withCache(boolean c = true) {
        cached = c
        return this
    }

    String getLabel() {
        Futils.removeExtention(name)
    }

    int getSize() {
        items.size()
    }

    /**
     * clear cached properties of cached proteins
     */
    void clearSecondaryCaches() {
        items.each {
            if (it.cachedPair!=null) {
                it.cachedPair.prediction.protein.clearSecondaryData()
                it.cachedPair.liganatedProtein.clearSecondaryData()
            }
        }
    }

    /**
     * clear cached structures
     */
    void clearPrimaryCaches() {
        items.each { it.cachedPair = null }
    }

    boolean checkFilesExist() {
        boolean ok = true
        items.each {
            if (header.contains(COLUMN_PREDICTION)) {
                if (!Futils.exists(it.pocketPredictionFile)) {
                    log.error "prediction file doesn't exist: $it.pocketPredictionFile"
                    ok = false
                }
            }
            if (header.contains(COLUMN_PROTEIN)) {
                if (!Futils.exists(it.proteinFile)) {
                    log.error "protein file doesn't exist: $it.proteinFile"
                    ok = false
                }
            }

        }
        return ok
    }

    /**
     * Process all dataset items with provided processor.
     */
    Result processItems(final Processor processor) {
        return processItems(params.parallel, processor)
    }

    /**
     * Process all dataset items with provided processor.
     */
    Result processItems(boolean parallel, final Processor processor) {

        log.trace "ITEMS: " + StrUtils.toStr(items)

        Result result = new Result()

        if (parallel) {
            int nt = ThreadPoolFactory.pool.poolSize
            log.info "processing dataset [$name] using $nt threads"

            AtomicInteger counter = new AtomicInteger(1);

            GParsPool.withExistingPool(ThreadPoolFactory.pool) {
                items.eachParallel { Item item ->
                    int num = counter.getAndAdd(1)
                    processItem(item, num, processor, result)
                }
            }

        } else {
            log.info "processing dataset [$name] using 1 thread"

            int counter = 1
            items.each { Item item ->
                processItem(item, counter++, processor, result)
            }
        }

        return result
    }

    private void processItem(Item item, int num, Processor processor, Result result) {

        String msg = "processing [$item.label] ($num/$size)"
        log.info (
           "\n------------------------------------------------------------------------------------------------------------------------"
         + "\n$msg"
         + "\n------------------------------------------------------------------------------------------------------------------------"
        )
        System.out.println(msg)

        try {

            processor.processItem(item)

        } catch (Exception e) {
            String emsg = "error processing dataset item [$item.label]"
            log.error(emsg, e)
            result.errorItems.add(item)

            if (params.fail_fast) {
                throw new PrankException(emsg, e)
            }
        }

    }

    private PredictionLoader getLoader(Item item) {
        return getLoader(attributes.get(PARAM_PREDICTION_METHOD), item)
    }

    private PredictionLoader getLoader(String method, Item item) {
        PredictionLoader res
        switch (method) {
            case "fpocket":
                res = new FPockeLoader()
                break
            case "concavity":
                res = new ConcavityLoader()
                break
            case "sitehound":
                res = new SiteHoundLoader()
                break
            default:
                res = new FPockeLoader()
                //throw new Exception("invalid method: $method")
        }

        if (res!=null) {
            res.loaderParams.ligandsSeparatedByTER = (attributes.get(PARAM_LIGANDS_SEPARATED_BY_TER) == "true")  // for bench11 dataset
            res.loaderParams.relevantLigandsDefined = hasLigandCodes()
            res.loaderParams.relevantLigandNames = item.getLigandCodes()
        }

        return res
    }

    public Dataset randomSubset(int subsetSize, long seed) {
        if (subsetSize >= this.size) {
            return this
        }

        List<Item> shuffledItems = new ArrayList<>(items)
        Collections.shuffle(shuffledItems, new Random(seed))

        return createSubset( shuffledItems.subList(0, subsetSize), this.name + " (random subset of size $subsetSize)" )
    }

    List<Fold> sampleFolds(int k, long seed) {
        if (size < k)
            throw new PrankException("There is less dataset items than folds! ($k < $size)")

        List<Item> shuffledItems = new ArrayList<>(items)
        Collections.shuffle(shuffledItems, new Random(seed))

        // split to subsets
        List<Item>[] subsets = new List<Item>[k]
        for (int i=0; i<k; i++)
            subsets[i] = new ArrayList<Item>()
        for (int i=0; i!=items.size(); ++i)
            subsets[i % k].add( shuffledItems[i] )

        // create folds
        List<Fold> folds = new ArrayList<>(k)
        for (int i=0; i!=k; ++i) {
            Dataset evalset = createSubset(subsets[i], "${name}_fold.${k}.${i}_eval")
            Dataset trainset = createSubset( shuffledItems - subsets[i], "${name}_fold.${k}.${i}_train" )
            folds.add(new Fold(i+1, trainset, evalset))
        }

        return folds
    }

    /**
     * create dataset view from subset of items
     */
    private Dataset createSubset(List<Item> items, String name) {
        Dataset res = new Dataset(name, this.dir)
        res.items = items
        res.attributes = this.attributes
        res.cached = this.cached
        res.header = this.header

        return res
    }

    /**
     *
     * @param pdbFile
     * @param itemContext allows to specify additional columns May be null.
     * @return
     */
    public static Dataset createSingleFileDataset(String pdbFile, ProcessedItemContext itemContext) {
        Dataset ds = new Dataset(Futils.shortName(pdbFile), Futils.dir(pdbFile))

        Map<String, String> columnValues = new HashMap<>()
        //columnValues.put(COLUMN_PROTEIN, pdbFile)
        //columnValues.put(COLUMN_PREDICTION, pdbFile)

        if (itemContext!=null) {
            columnValues.putAll(itemContext.datsetColumnValues)
        }
        
        ds.items.add(ds.createNewItem(pdbFile, null, columnValues))

        return ds
    }

    public static Dataset createSingleFileDataset(String pdbFile) {
        createSingleFileDataset(pdbFile, null)
    }

    public boolean hasLigandCodes() {
        return header.contains(COLUMN_LIGAND_CODES)
    }

    Dataset(String name, String dir) {
        this.name = name
        this.dir = dir
    }

     /**
     * file format:
     * <pre>
     * {@code
     *
     * # comment .. global dataset attributes start line with "PARAM."
     * PARAM.METHOD=fpocket/concavity
     *
     * # uncomment if ligands are separated with TER record:
     *  PARAM.LIGANDS_SEPARATED_BY_TER=true
     *
     * # relative paths
     * xxx/pocket-prediction-output1.pdb  liganateded-protein-file1.pdb
     * xxx/pocket-prediction-output1.pdb  liganateded-protein-file2.pdb
     *
     * }
     * </pre>
     * @param fname
     * @return
     */
    static Dataset loadFromFile(String fname) {
        File file = new File(fname)

        if (!file.exists()) {
            throw new PrankException("cannot find dataset file [$file.name]")
        }

        log.info "loading dataset [$file.absolutePath]"

        Dataset dataset = new Dataset(file.name, file.parent)

        for (String line in file.readLines()) {
            line = line.trim()
            if (line.startsWith("#") || line.isEmpty()) {
                // ignore comments and empty lines
            } else if (line.startsWith("PARAM.")) {
                String paramName = line.substring(line.indexOf('.') + 1, line.indexOf('=')).trim()
                String paramValue = line.substring(line.indexOf('=') + 1).trim()
                dataset.attributes.put(paramName, paramValue)
            } else if (line.startsWith("HEADER:")) {
                dataset.header = parseHeader(line)
            } else {
                dataset.items.add(dataset.parseItem(line))
            }
        }
        
        log.debug("dataset header: {}", dataset.header)

        if (!dataset.checkFilesExist()) {
            throw new PrankException("dataset contains invalid files")
        }

        return dataset
    }

    private Item parseItem(String line) {
        List<String> cols = SPLITTER.splitToList(line)
        Map<String, String> colValues = new HashMap<>()
        header.eachWithIndex { String col, int i ->
            colValues.put(col, cols[i])
        }
        return createNewItem(colValues)
    }

    static List<String> parseHeader(String line) {
        SPLITTER.splitToList(line).tail()
    }

}
