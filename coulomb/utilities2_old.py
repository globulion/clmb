# ------------------------------------------ #
#                UTILITIES                   #
# ------------------------------------------ #

from numpy import zeros, array, outer,\
                  float64

__all__=[]

def array_outer_product( A, B, result=None ):
    """ Compute the outer-product in the terminal two dimensions of the given arrays.
    If the result array is provided, the results are written into it.
    """
    #print A.shape[:-1],B.shape[:-1]
    assert( A.shape[:-1] == B.shape[:-1] )
    if result is None:
        result=zeros( A.shape+B.shape[-1:], dtype=A.dtype )
    if A.ndim==1:
        result[:,:]=outer( A, B )
    else:
        for idx in xrange( A.shape[0] ):
            array_outer_product( A[idx,...], B[idx,...], result[idx,...] )
        
    return result

def array_outer_product_1_2( A, B, result=None):
    """ Compute outer-product of two tensors, A(...1) and B(...2)
    according to the formula:
    Result_...ijk = A_...i * B_...jk
    Note: tensor A has its last  i axe  of dimension  3  !
    Note: tensor B has its last jk axes of dimensions 3x3!
    """ 
    if len(B.shape)>2:
       f = B.shape[0]
       result = array_outer_product( A,B.reshape(f,9) ).reshape(f,3,3,3)
    else:
       result = array_outer_product( A,B.reshape(9) ).reshape(3,3,3) 
    return result

def array_outer_product_2_1( B, A, result=None):
    """ Compute outer-product of two tensors, B(...2) and A(...1)
    according to the formula:
    Result_...ijk = B_...ij * A_...k
    Note: tensor B has its last ij axes of dimensions 3x3!
    Note: tensor A has its last  k axe  of dimension  3  !
    """ 
    if len(B.shape)>2:
       f = B.shape[0]
       result = array_outer_product( B.reshape(f,9),A ).reshape(f,3,3,3)
    else:
       result = array_outer_product( B.reshape(9),A ).reshape(3,3,3) 
    return result
