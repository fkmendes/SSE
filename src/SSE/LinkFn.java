package SSE;

/*
 * This interface returns f(x), where f(x) = y can be a macroevolutionary parameter
 * (e.g., lambda, mu) or initial values of D (from esDs)
 */
public interface LinkFn {

    boolean refreshParams();

    // x = qu trait, y = macroevol parameter or D (in esDs)
    default double[] getY(double[] x, double[] y, boolean ignoreRefresh) {
        return y;
    }

    // x = qu trait, y = macroevol parameter or D (in esDs)
    default double[] getY(double[] x, double[] y, int[] nLeftNRightFlanksHi, String spName, boolean ignoreRefresh) {
        return y;
    }

    String getLinkFnName();

}
