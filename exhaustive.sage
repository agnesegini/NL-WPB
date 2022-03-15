"""
This file contains the code for exhaustively computing the distribution of the distance.
"""
#load basic functions and packages
load("code_fun.sage")

#given as input two integers n and v r computes the spherically punctured Reed Muller code of order 1 of lenght k=(n choose v), and computes exhaustively the distribution of the distance among from such code among all the vectors of lenght k and hw k//2

def test_exhaustive(n,v,verbose=False):
  t0=time.time()
  print "gen code"
  M,C,k,A=gen_code_per(n,v)
  print "gen words", k,k//2 
#  li=all_fix_hw(k,k//2)
#  print "list size: ", len(li)
  O=vector(ZZ,k//2)
  MAX=0
  it=ithp(k,k//2)
  for rp in it:
      r=vector(GF(2),k)
      for i in rp: r[i]=1
      D=dist_all(C,vector(r))
      h=min(D)
      if verbose: print r,h#,D
      del D,r
      O[h]+=1
      if h>MAX: MAX=h
  print '--->MAX distance :', MAX 
  print "\nrunning time: %.2f" %(time.time()-t0)
  pl=plot_bar(O,name,per=False,color=orange):
  return O,MAX,pl
   
   
#loop for the parallelisation of the exhaustive search, see definition of test_exhaustive_para 
def loop_ex(k,n,C,dist_,rp,verbose=False):
  r=vector(GF(2),k)
  for i in rp: r[i]=1
  D=dist_(C,r)
  h=min(D)
  del D   
#  if h>26 or verbose: print h, r
  return h
  
#given as input two integers n and v, computes the spherically punctured Reed Muller code of order 1 of length k=(n choose v), and computes exhaustively the distribution of the distance among from such code among all the vectors of lenght k and hw k//2. The seach is performed in parallel over 8 workers.   The output is the a vector O such that O[i] is the number of string at distance i, a  plot pl of suc distribution and MAX the max distance found. 
def test_exhaustive_para(n,v,verbose=False):
    t0=time.time()
    print "gen code"
    M,C,k,A=gen_code_per(n,v)
    b=binomial(k,k//2)
    print "gen words", k,k//2,  binomial(k,k//2)
    O=vector(ZZ,k//2)
    MAX=0
    it=ithp(k,k//2)
    dist_=dist_all
    pool = Pool(8) # start 8 worker processes
    fo=partial(loop_ex,k,n,C,dist_)
    H=pool.map(fo,it)
    assert len(H)==b
    pool.close()
    pool.join()
    for h in H:
      if h>MAX: MAX=h
      O[h]+=1
    print '--->MAX distance :', MAX 
    print "\nrunning time: %.2f" %(time.time()-t0)
    pl=plot_bar(O,name,per=False,color=orange):
    return O,MAX,pl
    
    

