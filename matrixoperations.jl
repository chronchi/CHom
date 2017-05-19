#row exchange between i and j rows

function rowExchange(B::Matrix, i, j)
  B = copy(B)
  B[i,:],B[j,:] = B[j,:], B[i,:]
  return convert(Array{Int},B)
end

#multiply the i-th row by -1

function rowMultiply(B::Matrix, i)
  B = copy(B)
  B[i,:] = -B[i,:]
  return convert(Array{Int},B)
end

#substitute the i-th row by the i-th row plus q times the j-th row

function rowAdd(B::Matrix, i::Int, j, q)
  if typeof(q) != Int
    error("q must be an integer")
    return
  end
  n = size(B,2)
  B = copy(B)
  B[i,1:n] = B[i,1:n] + q*B[j,1:n]
  return convert(Array{Int},B)
end

#exchange columns i and j

function columnExchange(B::Matrix, i, j)
  B = copy(B)
  B[:,i], B[:,j] = B[:, j], B[:, i]
  return convert(Array{Int},B)
end

#multiply column i by -1

function columnMultiply(B::Matrix, i)
  B = copy(B)
  B[:,i] = -B[:,i]
  return convert(Array{Int},B)
end

#substitute the i-th column by the i-th column plus q times the j-th column

function columnAdd(B::Matrix, i::Int, j, q)
  m = size(B,1)
  B = copy(B)
  B[1:m, j] = B[1:m, j] + q*B[1:m, i]
  return convert(Array{Int},B)
end

#Row exchange operation keeping track of bases

function rowExchangeOperation(B,Q,Qbar,i,j)
  B = rowExchange(B,i,j)
  Qbar = rowExchange(Qbar,i,j)
  Q = columnExchange(Q,i,j)
  return B,Q,Qbar
end

#Multiply a row by -1 keeping track of bases

function rowMultiplyOperation(B,Q,Qbar,i)
  B = rowMultiply(B,i)
  Qbar = rowMultiply(Qbar,i)
  Q = columnMultiply(Q,i)
  return B,Q,Qbar
end

#Add a multiple of a row keeping track of bases

function rowAddOperation(B,Q,Qbar,i,j,q)
  B = rowAdd(B,i,j,q)
  Qbar = rowAdd(Qbar, i, j, q)
  Q = columnAdd(Q,i,j,-q)
  return B,Q,Qbar
end

#Column exchange operation keeping track of bases

function columnExchangeOperation(B,R,Rbar,i,j)
  B = columnExchange(B,i,j)
  Rbar = columnExchange(Rbar,i,j)
  R = rowExchange(R,i,j)
  return B,R,Rbar
end

#Multiply a row by -1 keeping track of bases

function columnMultiplyOperation(B,R,Rbar,j)
  B = columnMultiply(B,j)
  Rbar = columnMultiply(Rbar,j)
  R = rowMultiply(R,j)
  return B,R,Rbar
end

#Add a multiple of a column keeping track of bases

function columnAddOperation(B,R,Rbar,i,j,q)
  B = columnAdd(B,i,j,q)
  Rbar = columnAdd(Rbar,i,j,q)
  R = rowAdd(R,i,j,-q)
  return B,R,Rbar
end

#Algorithm 3.29 Partial row reduction

function partRowReduce(B,Q,Qbar,k,l)
  for i = (k+1):size(B,1)
    q = convert(Int64, floor(B[i,l]/B[k,l]))
    B,Q,Qbar = rowAddOperation(B,Q,Qbar,i,k,-q)
  end
  return B,Q,Qbar
end

#partial column reduction

function partColumnReduce(B,R,Rbar,k,l)
  for i = (k+1):size(B,2)
    q = convert(Int64, floor(B[l,i]/B[l,k]))
    B, Q, Qbar = columnAddOperation(B, Q, Qbar, i, k, -q)
  end
  return B, Q, Qbar
end

#Algorithm 3.32 Smallest nonzero entry

function smallestNonzero(v,k)
  c = length(v)
  alpha = minimum([abs(v[i]) for i=k:c if v[i] != 0])
  i0 = minimum([i for i=k:c if abs(v[i]) == alpha])
  return alpha, i0
end

#Algorithm 3.33 Row preparation

function rowPrepare(B,Q,Qbar,k,l)
  alpha, i = smallestNonzero(B[:,l], k)
  B,Q,Qbar = rowExchangeOperation(B,Q,Qbar,k,i)
  return B,Q,Qbar
end

#Algorithm 3.36 Row reduction

function rowReduce(B,Q,Qbar,k,l)
  m = size(B,1)
  while B[(k+1):m,l] != zeros(Int, m-k)
    B,Q,Qbar = rowPrepare(B,Q,Qbar,k,l)
    B,Q,Qbar = partRowReduce(B,Q,Qbar,k,l)
  end
  return B,Q,Qbar
end

#Algorithm 3.39 Row echelon

function rowEchelon(A)
  m,n = size(A)
  B = copy(A)
  Q = eye(Int, m)
  Qbar = eye(Int, m)
  l = 1
  k = 0
  while k < m
    while l <= n && B[(k+1):m,l] == zeros(Int, m-k)
      l += 1
      println(l)
    end
    if l == n+1
      break
    end
    k += 1
    B,Q,Qbar = rowReduce(B,Q,Qbar,k,l)
    println("finished row reducing")
  end
  return B,Q,Qbar,k
end

#Algorithm 3.42 Kernel-image algorithm

function kernelImage(B)
  m,n = size(B)
  BT = B'
  B,P,Pbar,k = rowEchelon(BT)
  BT = B'
  PT = P'
  return PT[:,(k+1):n], BT[:,1:k]
end
