package mosse;

import beast.core.parameter.RealParameter;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertEquals;

public class MosseTipLikelihoodTest {

    private static double DELTA = 1e-10;

    @Test
    public void testMosseTipLikelihoodLargeIntervalIsOne() {
        double beta0 = 0.1;
        double beta1 = 0.2;
        double epsilon = 0.01;

        double[] traits = {1.0, 1.0};
        double a = -5.0;
        double b = 5.0;

        Double[] betasArray = {beta0, beta1};
        RealParameter betas = new RealParameter(betasArray);

        MosseTipLikelihood tipLikelihood = new MosseTipLikelihood();
        tipLikelihood.initByName(
                "beta", betas,
                "epsilon", Double.toString(epsilon)
        );
        tipLikelihood.initAndValidate();

        double prob = tipLikelihood.getTipLikelihood(a, b, traits);

        assertEquals(1.0, prob, DELTA);
    }
}