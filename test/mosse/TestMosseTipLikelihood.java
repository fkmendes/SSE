package mosse;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class TestMosseTipLikelihood {

    private static double DELTA = 1e-10;

    @Test
    public void testMosseTipLikelihoodLargeIntervalIsOne() {
        double beta0 = 0.1;
        double beta1 = 0.2;
        double epsilon = 0.01;

        double trait0 = 1.0;
        double trait1 = 1.0;
        double a = -5.0;
        double b = 5.0;

        MosseTipLikelihood tipLikelihood = new MosseTipLikelihood();
        tipLikelihood.initByName(
                "beta0", Double.toString(beta0),
                "beta1", Double.toString(beta1),
                "epsilon", Double.toString(epsilon)
        );
        tipLikelihood.initAndValidate();

        double prob = tipLikelihood.getTipLikelihood(a, b, trait0, trait1);

        assertEquals(1.0, prob, DELTA);
    }
}