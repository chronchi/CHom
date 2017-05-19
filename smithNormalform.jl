include("rowEchelonform.jl")

#Algorithm 3.46 Minimal nonzero entry

#function minNonzero(B,k)
#  v = []
#  q = []
#  m,n = size(B)
#  for i = 1:m
#    if i < k
#      v = [v;0]
#      q = [q;0]
#    else
#      b = smallestNonzero(B[i,:],k)[1]
#      v = [v;b]
#      q = [q;b]
#    end
#  end
#  alpha, j = smallestNonzero(v,k)
#  return alpha, j, q[j]
#end

function minNonzero(B,k)
  if k == size(B,2)
    return B[k,k], k, k
  end
  C = abs.(B[k:end,k:end])
  d = zeros(Int, 0,2)
  for i = 1:(n-k+1)
    a,b = smallestNonzero(C[:,i],1)
    d = [d; a b]
  end
  e,f = smallestNonzero(d[:,1],1) #menor valor e coluna do menor valor
  g = d[f,2] #linha do menor valor
  return e, f+k-1, g+k-1
end

#Algorithm 3.47 Move minimal nonzero entry

function moveMinNonzero(B,Q,Qbar,R,Rbar,k)
  alpha, i, j = minNonzero(B,k)
  B,Q,Qbar = rowExchangeOperation(B,Q,Qbar,k,j)
  B,R,Rbar = columnExchangeOperation(B,R,Rbar,k,i)
  return B,Q,Qbar,R,Rbar
end

#Algorithm 3.48

function checkForDivisibility(B,k)
  m,n = size(B)
  for i = (k+1):m
    for j = (k+1):n
      q = floor(Int,B[i,j]/B[k,k])
      if q*B[k,k] != B[i,j]
        return false, i, j, q
      end
    end
  end
  return true,0,0,0
end

#Algorithm 3.49 which is wrong

function partSmithForm(B,Q,Qbar,R,Rbar,k)
  m,n = size(B)
  divisible = false
  while divisible == false
    B,Q,Qbar,R,Rbar = moveMinNonzero(B,Q,Qbar,R,Rbar,k)
    B,Q,Qbar = partRowReduce(B,Q,Qbar,k,k)
    if B[k+1:m,k] != zeros(Int, m-k)
      continue
    end
    divisible, i, j, q = checkForDivisibility(B,k)
    if divisible == true
      break
    end
    println("oi linha 75")
    B,R,Rbar = partColumnReduce(B,R,Rbar,k,k)
    if B[k,k+1:n] != zeros(Int, n-k)
      continue
    end
    println("linha 80")
    divisible, i, j, q = checkForDivisibility(B,k)
    println(divisible)
    if divisible == false
      B,Q,Qbar = rowAddOperation(B,Q,Qbar,i,k,q)
      B,R,Rbar = columnAddOperation(B,R,Rbar,k,j,-q)
    end
  end
  return B,Q,Qbar,R,Rbar
end

#Algorithm 3.51 Smith algorithm

function smithForm(B)
  m,n = size(B)
  Q = Qbar = eye(Int, m)
  R = Rbar = eye(Int, n)
  s = t = 0
  while B[t+1:m,t+1:n] != 0
    t += 1
    B,Q,Qbar,R,Rbar = partSmithForm(B,Q,Qbar,R,Rbar,t)
    if B[t,t] < 0
      B,Q,Qbar = rowMultiplyOperation(B,Q,Qbar,t)
    end
    if B[t,t] == 1
      s += 1
    end
  end
  return B,Q,Qbar,R,Rbar,s,t
end


#B = A;
#Q = eye(Int,3);
#Qbar = eye(Int,3);
#R = eye(Int,3);
#Rbar = eye(Int,3);
