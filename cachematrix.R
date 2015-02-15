## Defines a special matrix object (makeCacheMatrix) that can cache its 
## inverse, and a function 'cacheSolve' that computes and sets the inverse. 

## Creates an object that stores a matrix and caches its inverse. 
##      Args: 
##          A 'matrix'. 
##      Returns: 
##          The 'list' of functions: 
##          - get(), returns the stored matrix; 
##          - set(), redefines the stored matrix; 
##          - isSquared(), whether it is a squared matrix; 
##          - setInv(), caches the matrix's inverse; 
##          - getInv(), returns the matrix's inverse.
makeCacheMatrix <- function(mat = matrix(c(1,0,0,1),2,2)){
    matInv <- NULL; # matrix's inverse
    squared <- ncol(mat)==nrow(mat); # is squared matrix?
    set <- function(matp) {
        mat <<- matp;
        matInv <<- NULL;
        squared <<- ncol(matp)==nrow(matp);
    };
    get <- function() return(mat);
    isSquared <- function() return(squared);
    setInv <- function(matpInv) matInv <<- matpInv;
    getInv <- function() return(matInv);
    return(list(set = set, get = get,
         setInv = setInv,
         getInv = getInv, isSquared = isSquared));
}

## Computes and caches the inverse of a matrix stored in a 'makeCacheMatrix' object. 
## (The computation is skipped if the inverse matrix is already cached in the object.) 
## If the matrix is non-squared, the Moore-Penrose pseudoinverse is considered. 
##      Args: 
##          A 'makeCacheMatrix' object. 
##      Returns: 
##          A 'matrix' that is the inverse of the matrix stored in the argument. 
cacheSolve <- function(mat = makeCacheMatrix()){
    if (is.null(mat$getInv())){ 
        if (mat$isSquared())
            matInv<-solve(mat$get())
        else { 
            message("Non-squared matrix: calculating the Moore-Penrose pseudoinverse");
            # calculate pseudoinverse by the SVD method
            duv<-svd(mat$get());
            matInv<-duv$v%*%Conj(diag(1/duv$d,length(duv$d)))%*%Conj(t(duv$u));
        };
        mat$setInv(matInv);
    };
    return(mat$getInv());
}