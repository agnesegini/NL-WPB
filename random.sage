"""
This file contains the code for randomly computing the distribution of the distance.
"""

#load basic functions and packages
load("code_fun.sage")

#given as input two integers n and v returns a vector w of length n and Hamming weight v. Remark: using the function randint_par as defined in code_fun.sage is essential for using this function inside a pool.
def gen_vec(n, v):
  l=0
  L=[]
  w=vector(GF(2),n)
  while l<v:
      p=randint_par(n)
      if p not in L: 
        L+=[p]
        l+=1
        w[p]=1     
  return w


   
#loop for the parallelisation of the random search, see definition of test_rand_para 
def loop_rand(k,n,v,C,dist_,I,BOUND=matrix_max):
  w=gen_vec(k,k//2)
  D=dist_(C,w)
  h=min(D)
  del D   
  if h>=BOUND[n,v]: print '\n', h, w
  return h,w


#given as input two integers n and v, computes C the spherically punctured Reed Muller code of order 1 of lenght k=(n choose v); and computes sequence, of lenght coeff_par, of integers, such that each one of them is the distance between C and a random sampled vector of lenght k and hw k//2 (this computation is performed in parallel over 8 workers); then the local max (resp. min) distance is compared with the max (resp. min) computed and in case this is higher (resp. lower) updated. n_sample//coeff_par are analised and compared. The output are the max and min distance found over (n_sample//coeff_par)*coeff_par samples.

def test_rand_para(n,v,n_sample=1000,coeff_par=24,verbose=False):
    t0=time.time()
    print "gen code"
    M,C,k,A=gen_code_per(n,v)
    b=binomial(k,k//2)
    print "gen words", k,k//2, b
    MIN=b
    MAX=0
    treasure=[]
    i=0
    dist_=dist_all
    n_iter=int(n_sample/coeff_par)
    while i<n_iter:
        pool = Pool() #on my laptop is 8
        fo=partial(loop_rand,k,n,v,C,dist_)
        LV=pool.map(fo,xrange(coeff_par))
        pool.close()
        pool.join()
        H=[c[0] for c in LV]
        MAX_loc=max(H)
        MIN_loc=min(H)
        #if MAX_loc>26: print MAX_loc,LV[H.index(MAX_loc)][1]
        del LV,H
        i+=1
        if MAX_loc>MAX: MAX=MAX_loc
        if MIN_loc<MIN: MIN=MIN_loc
        print (MAX_loc,MIN_loc),
        
    print i
    print '--->MIN distance :', MIN 
    print '--->MAX distance :', MAX 
    print "\nrunning time: %.2f" %(time.time()-t0)
    return MAX,MIN
    
    
   
#loop for the parallelisation of the random distribution estimation, see definition of test_rand_para 
def loop_dist_rand(k,n,v,C,dist_,I,BOUND=matrix_max):
  w=gen_vec(k,k//2)
  D=dist_(C,w)
  h=min(D)
  del D   
  if h>=BOUND[n,v]: print '\n', h, w
  return h,w


##given as input two integers n and v, computes C the spherically punctured Reed Muller code of order 1 of lenght k=(n choose v); and computes sequence, of lenght coeff_par, of integers, such that each one of them is the distance between C and a random sampled vector of lenght k and hw k//2 (this computation is performed in parallel over 8 workers); then the local max (resp. min) distance is compared with the max (resp. min) computed and in case this is higher (resp. lower) updated. n_sample//coeff_par are analised and compared. The output is an estimation of the distribution D_{v,n}.
def distribution_rand_para(n,v,n_sample=1000,coeff_par=24,verbose=False,outfile=None,BOUND=matrix_max):
    t0=time.time()
    print "gen code"
    M,C,k,A=gen_code_per(n,v)
    b=binomial(k,k//2)
    print "gen words", k,k//2, b
    MIN=k
    MAX=0
    O=vector(ZZ,k//2)
    treasure=[]
    i=0
    dist_=dist_all
    n_iter=int(n_sample/coeff_par)
    if outfile!=None:
      #f = open('/home/agnese/WPB/nlWPB/code'+outfile, 'w')
      #f = open('simul/'+outfile, 'w')
      f = open(outfile, 'w')
      print 'Output file:', f
      print outfile
      f.write("---Random sample distribution---\n")
      f.write("n: "+str(n)+",k :"+str(v)+",n_sample: "+str(n_sample)+'\n')     
    while i<n_iter:
        MAX_loc=0
        MIN_loc=k
        pool = Pool(8)
        fo=partial(loop_rand,k,n,v,C,dist_)
        LV=pool.map(fo,xrange(coeff_par))
        pool.close()
        pool.join()
        for c in LV: 
              h=c[0]
              O[h]+=1
              if h>MAX_loc: MAX_loc=h
              if MIN_loc>h: MIN_loc=h
              if (outfile!=None) and (c[0]>=BOUND[n,v]): f.write('\n'+str(h)+'--'+str(c[1]))
        #if MAX_loc>26: print MAX_loc,LV[H.index(MAX_loc)][1]
        del LV
        i+=1
        if MAX_loc>MAX: MAX=MAX_loc
        if MIN_loc<MIN: MIN=MIN_loc
        print (MAX_loc,MIN_loc),
    t1= time.time()   
    print i
    print '--->MIN distance :', MIN 
    print '--->MAX distance :', MAX 
    print "\nrunning time: %.2f" %(t1-t0)
    if outfile!=None:
      f.write("\n\nO: "+str(O)+"\nMAX :"+str(MAX)+",MIN: "+str(MIN))
      f.write("\nrunning time: %.2f sec" %(t1-t0))
      f.close()
    
    return O,MAX
    
def testrandom(n,v,nsa):
  name='simul/'+str(n)+'_'+str(v)+'**'+now_str()
  O,MAX=distribution_rand_para(n,v,nsa,outfile=name+'.txt')
  plot_bar(O,name,blue)
