package cz.siret.prank.domain

import cz.siret.prank.domain.labeling.LabeledPoint
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Pocket prediction result for single protein.
 */
@Slf4j
@CompileStatic
class Prediction {

    Protein protein

    /**
     * pockets predicted by P2RANK or other prediction method
     */
    List<Pocket> pockets

    /**
     * reordered pockets (relevant only when doing rescoring with old PRANK algorihhm)
     */
    List<Pocket> reorderedPockets

    /**
     *  SAS points with ligandability score for prediction and visualization.
     */
    List<LabeledPoint> labeledPoints = null


    Prediction(Protein protein, List<? extends Pocket> pockets) {
        this.protein = protein
        this.pockets = (List<Pocket>) pockets
    }

    int getPocketCount() {
        return pockets.size()
    }

}
