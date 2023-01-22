package mosse;

import beast.core.Function;
import beast.core.Input;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

/**
 * @author Kylie Chen
 * setup custom Q matrix for testing Mosse Tree Likelihood only
 */
public class CustomSubstitutionModel extends GeneralSubstitutionModel {


    final public Input<Function> customRatesInput =
            new Input<>("customRates", "Rate parameter which defines the transition rate matrix exactly. ", Input.Validate.REQUIRED);


    @Override
    public void initAndValidate() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        super.initAndValidate();
    }
    @Override
    protected void setupRateMatrixUnnormalized() {
        double[] customRatesValues = new double[nrOfStates * nrOfStates];
        Function customRates = this.customRatesInput.get();
        for (int i = 0; i < customRates.getDimension(); i++) {
            customRatesValues[i] = customRates.getArrayValue(i);
        }
        int count = 0;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = customRatesValues[count];
                count++;
            }
        }
        System.out.println();
    }
}
