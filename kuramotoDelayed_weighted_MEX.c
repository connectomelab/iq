#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

int mainf(float* theta, unsigned tmax,
          unsigned numItersPerSec, unsigned char numNodes,
          float k, unsigned short* m_irs, unsigned short* m_jcs, 
          float* mFLOAT,
          unsigned short* dINT, float* oscPI, unsigned tBuffer,
          float* pastStates);

int mainf(float* theta, unsigned tmax,
          unsigned numItersPerSec, unsigned char numNodes,
          float k, unsigned short* m_irs, unsigned short* m_jcs,
          float* mFLOAT,
          unsigned short* dINT, float* oscPI, unsigned tBuffer,
          float* pastStates){

    float dt = 1/(float)numItersPerSec;
    unsigned tmaxmax = tmax+tBuffer;
    unsigned t;
    float noise;
    float interactionTerm;
    for (t=tBuffer; t<tmaxmax; t++) {
        //printf("tBuffer:%d, time:%d/%d \n", tBuffer, t, tmax);
        unsigned curCol = t * numNodes;
        unsigned i;
        for (i=numNodes; i>0; i--){
           // printf(" %d/%d \n", i, numNodes);
            // num vals encountered so far
            unsigned short nnzSoFar = m_jcs[i-1];
            // num entries in column
            unsigned short numConns = m_jcs[i]-nnzSoFar;
            // pointer to connected node row-index (first in column)
            unsigned short* startOfCol_firstNZRowIdx = m_irs+nnzSoFar;
            // current theta value for this node
            float ti = theta[i-1];
            
            // for each connection into this node
            interactionTerm = 0;
            unsigned short connectedToNodeIdx;
            unsigned delayedColIdx; // leave unsigned
            float connectionWeight;
            unsigned x;
            for (x=0;x<numConns;x++){
                // column index into the vector of previous thetas
                delayedColIdx = (curCol-dINT[nnzSoFar+x]);

                // index of connected node
                connectedToNodeIdx = *(startOfCol_firstNZRowIdx+x);
                
                // log of number of streamlines (weight)
                connectionWeight = mFLOAT[nnzSoFar+x];

                // weighted difference
                interactionTerm += connectionWeight *
                    sin(pastStates[delayedColIdx + connectedToNodeIdx] - ti);
            }

            // NOTE: uncomment line to take AVERAGE input across all incoming connections
            // - watch for division by zero
            // interactionTerm = interactionTerm / numConns
        
            // Update current theta
            theta[i-1] = theta[i-1] + (dt*(k*(interactionTerm)+oscPI[i-1]));

            // add to previous thetas (for delayed theta values)
            pastStates[curCol+(i-1)] = theta[i-1];
        }
//         if (i<0.5){
//             printf("------t=%d\n", t);
//             printf(" osc: %f\n", oscPI_dt);
//             printf(" k*interaction: %f\n", k*interactionTerm*dt);
//         }
    }
    //printf("\n");
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

    float *initial, k, tsec, *oscPI;
    // Sparse has to be double.
    double *matrix, *delays;
    mwIndex *m_irs,*m_jcs;
    unsigned tBuffer, tmax, numItersPerSec;
    unsigned char numNodes;
    mwSize numNZ;
    
    /////////////////
    //// Validate
    /////////////////
    if(!mxIsSingle(prhs[0]) ||
       mxIsScalar(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:myfcn:nrhs",
                          "theta init must be single, non-scalar.");
    }
    if(!mxIsDouble(prhs[1]) ||
       !mxIsScalar(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:myfcn:nrhs",
                          "time (sec) must be double, scalar.");
    }
    if(!mxIsDouble(prhs[2]) ||
       !mxIsScalar(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:myfcn:nrhs",
                          "iters per sec must be double, scalar.");
    }
    if(!mxIsSparse(prhs[3]) ||
       !mxIsDouble(prhs[3])) {
        mexErrMsgIdAndTxt("MyToolbox:myfcn:nrhs",
                          "conn. matrix must be sparse double.");
    }
    if(!mxIsSparse(prhs[4]) ||
       !mxIsDouble(prhs[4])) {
        mexErrMsgIdAndTxt("MyToolbox:myfcn:nrhs",
                          "delays must be sparse double.");
    }
    if(!mxIsSingle(prhs[5]) ||
       !mxIsScalar(prhs[5])) {
        mexErrMsgIdAndTxt("MyToolbox:myfcn:nrhs",
                          "global weight must be single, scalar.");
    }
    if(!mxIsSingle(prhs[6]) ||
       mxIsSparse(prhs[6]) ||
       mxIsScalar(prhs[6])) {
        mexErrMsgIdAndTxt("MyToolbox:myfcn:nrhs",
                          "osc values must be single, non-sparse, non-scalar.");
    }
    if(!mxIsSingle(prhs[7]) ||
       !mxIsScalar(prhs[7])) {
        mexErrMsgIdAndTxt("MyToolbox:myfcn:nrhs",
                          "max delay must be single, scalar.");
    }

    /////////////////
    //// Assign
    /////////////////
    // starting values for each node (theta)
    initial = mxGetSingles(prhs[0]);
    // duration in seconds of simulation
    tsec = (float)mxGetScalar(prhs[1]);
    // number of iterations per second
    numItersPerSec = (unsigned)(mxGetScalar(prhs[2])+0.5);

    /// Connection matrix; extract useful bits
    // number of nodes (e.g. 82)
    numNodes = (unsigned char)mxGetM(prhs[3]);
    // number of connections
    numNZ = mxGetNzmax(prhs[3]);
    // connection matrix (log of # of streamlines)
    matrix = mxGetDoubles(prhs[3]);
    // irs: contains nnz integers (rows [offset by 1] where non-zero elements are to be found)
    m_irs = mxGetIr(prhs[3]); 
    // m_jcs[j]: the index into irs of the first non-zero value for that jth column
    // The last element of the jc array, jc[n], is equal to nnz: number of nonzero elements in sparse matrix.
    m_jcs = mxGetJc(prhs[3]);
    
    // delays (time-steps when looking at past values)
    delays = mxGetDoubles(prhs[4]);
    // global weight multipled against each connection
    k = (float)mxGetScalar(prhs[5]);
    // each node will have a slightly different oscillation freq (noise applied earlier)
    oscPI = mxGetSingles(prhs[6]);
    // maximum delay
    tBuffer = (unsigned) (mxGetScalar(prhs[7])+0.5);

    //// Output (PAST STATES)
    tmax = (unsigned) (tsec*numItersPerSec + 0.5);
    
    float* theta = (float*)malloc(numNodes*sizeof(float));
    memcpy(theta, initial, sizeof(float)*numNodes);
    
    // we want pr to hold floats (computational savings)
    unsigned short* dINT = (unsigned short*)malloc(numNZ*sizeof(unsigned short));
    unsigned short* irsINT = (unsigned short*)malloc(numNZ*sizeof(unsigned short));
    float* mFLOAT = (float*)malloc(numNZ*sizeof(float));
    int i;
    for (i=0;i<numNZ;i++){
        dINT[i] = (unsigned short)(delays[i]+0.5) * numNodes; // num nodes = jump when looking backwards
        irsINT[i] = (unsigned short)m_irs[i];
        mFLOAT[i] = (float)matrix[i];
        //printf("%f %u\n", delays[i], dINT[i]);
    }
    
    unsigned short* jcsINT = (unsigned short*)malloc((numNodes+1)*sizeof(unsigned short));
    for (i=0;i<numNodes+1;i++){
        jcsINT[i] = (unsigned short)m_jcs[i];
    }
    
    printf("MEX: numNodes:%d, tSize:%d, tsec:%f, numItersPerSec:%d, k:%f\n",
            numNodes, tmax+tBuffer, tsec, numItersPerSec, k);
    

    /////////////
    plhs[0] = mxCreateNumericMatrix(numNodes, (mwSize) (tmax+tBuffer), mxSINGLE_CLASS, mxREAL);
    float* pastStates = mxGetSingles(plhs[0]);
    mainf(theta, tmax, numItersPerSec, numNodes,
         k, irsINT, jcsINT, mFLOAT, dINT, oscPI, tBuffer, pastStates);

    free(dINT);
    free(theta);
    free(irsINT);
    free(jcsINT);
}
