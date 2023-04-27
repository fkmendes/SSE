package mosse;

import beast.core.Function;
import beast.core.Input;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

import java.lang.reflect.InvocationTargetException;

/**
 * @author Kylie Chen
 * setup custom Q matrix for testing Mosse Tree Likelihood only
 */
public class CustomSubstitutionModel extends GeneralSubstitutionModel {

    final public Input<Function> customRatesInput =
            new Input<>("customRates", "rate parameter which defines the transition rate matrix Q exactly", Input.Validate.REQUIRED);

    public CustomSubstitutionModel() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        frequencies = frequenciesInput.get();
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;

        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
                 | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
        rateMatrix = new double[nrOfStates][nrOfStates];
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
    }

    @Override
    protected void setupRelativeRates() {
        // use unnormalized matrix method instead
    }
}
