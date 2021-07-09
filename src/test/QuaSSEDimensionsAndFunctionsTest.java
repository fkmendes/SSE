package test;

import SSE.LogisticFunction;
import org.junit.Test;
import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

public class QuaSSEDimensionsAndFunctionsTest {

    final static double EPSILON = 1e-6;

    LogisticFunction lf;

    /*
     * Applies logistic function to many x values, checks return
     */
    @Test
    public void testLogistic() {

        double xInitial, x0, y1, y0, r, dx;
        double[] x = new double[3999];

        dx = 0.025;
        xInitial = 0.0;
        for (int i=0; i<4000; i++) {
            x[i] = xInitial;
            xInitial += dx;
        }


        // assertArrayEquals(expectedStartFy, Arrays.copyOfRange(fY, 0, 10), EPSILON);

    }
}
