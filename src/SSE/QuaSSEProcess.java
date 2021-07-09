package SSE;

import beast.core.Description;

@Description("Specifies a quantitative trait(s) state-dependent speciation and" +
        "extinction birth-death process.")
public interface QuaSSEProcess {

    /*
     *
     */
    void initializeTips();

    /*
     *
     */
    void pruneTree();

    /*
     *
     */
    void processBranch();

    /*
     *
     */
    void processInternalNode();

    /*
     *
     */
    void processRootNode();

    /*
     *
     */
    void doIntegrate();

    /*
     *
     */
    void propagateT();

    /*
     *
     */
    void propagateX();

    /*
     *
     */
    void convolve();
}
