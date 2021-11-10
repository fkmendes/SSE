package test;

import SSE.SSEUtils;
import org.junit.Test;

import javax.management.RuntimeErrorException;
import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

/**
 * @author Fabio K. Mendes
 */

/*
 * This class contains unit tests for
 */
public class OtherUtilsTest {

    final static double EPSILON = 1e-6;

    /*
     * Normalizes double array
     */
    @Test
    public void testNormalize() {
        double[] anArray = new double[] { 2.0, 3.0, 4.0 };
        double[] expected = new double[] { 0.222222, 0.333333, 0.444444 };
        SSEUtils.normalizeArray(anArray);

        assertArrayEquals(expected, anArray, EPSILON);
    }

    @Test
    public void testGetEveryOtherElement() {
        // input arrays
        double[] fromArray1 = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };
        double[] fromArray2 = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0 };

        // expected arrays
        /*
         * 0 = 1.0 (whatever was there), (skipped source 1.0), then ignore (2.0)
         * 1 = get (3.0), then ignore (4.0)
         * 2 = get (5.0), then ignore (6.0)
         * 3 = get (7.0), then ignore (8.0)
         * 4 = get (9.0)
         * 5... = 6.0... whatever was there
         */
        double[] toArrayExpected1 = new double[] {
                1.0, 3.0, 5.0, 7.0, 9.0, 6.0, 7.0, 8.0,
                9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };
        /*
         * 0 = 1.0 (whatever was there), (skipped source 1.0)
         * 1 = 2.0 (whatever was there), (skipped source 3.0)
         * 2 = get (5.0), then ignore (6.0)
         * 3 = get (7.0), then ignore (8.0)
         * 4 = get (9.0), then ignore (10.0)
         * 5 = get (11.0),
         * 6... = 7.0... whatever was there
         */
        double[] toArrayExpected2 = new double[] {
                1.0, 2.0, 5.0, 7.0, 9.0, 11.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0 };

        SSEUtils.everyOtherInPlace(fromArray1, 8, 1, 1, 2, 1.0);
        SSEUtils.everyOtherInPlace(fromArray2, 12, 2, 2, 2, 1.0);

        assertArrayEquals(toArrayExpected1, fromArray1, 0.0);
        assertArrayEquals(toArrayExpected2, fromArray2, 0.0);
    }

    /*
     * Will change this later because these checks will happen
     * in 'initAndValidate'. This verifies invalid inputs.
     */
    @Test(expected = RuntimeException.class)
    public void testCopyingTooManyElements() {

        // input arrays
        double[] inArray = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };

        // output arrays
        double[] outArray = new double[16];

        SSEUtils.everyOtherInPlace(inArray, 8, 2, 2, 2, 1.0);
    }

    /*
     * Test for SSE function that copies elements from high-res array of E's or D's
     * to low-res array. This function requires an array of indices to make the copy,
     * which is initialized and populated by the likelihood class upon initialization.
     */
    @Test
    public void testFromHiLoTransfer() {
        // input arrays
        double[] fromArray = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };
        double[] toArray = new double[4];
        int[] hiLoIdxs4Transfer = new int[] { 3, 7, 11, 15 };

        SSEUtils.hiToLoTransferInPlace(fromArray, toArray, hiLoIdxs4Transfer);

        double[] toArrayExpected = new double[] { 4.0, 8.0, 12.0, 16.0 };
        assertArrayEquals(toArrayExpected, toArray, 0.0);
    }
}
