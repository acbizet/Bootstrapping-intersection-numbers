"""Local exact-prefix height optimization on physical coefficient-torus angles.

This is a finite multistart diagnostic, not a certified global bound.  The
independent angles are theta_x,theta_y,u, all in radians.  No source saddle
condition or fixed-lambda relation is imposed while varying those angles.
"""
from pathlib import Path
import argparse,json,sys,time,re
import numpy as np
from scipy.optimize import brentq,minimize
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from product_engine import PrefixProduct
OUT=Path(__file__).resolve().parent
A=np.array([[1j,0,0],[-1j,3j,0],[0,0,1]],complex)

def wrap(x):return float(np.arctan2(np.sin(x),np.cos(x)))
def encode(S):return {'re':float(S.real),'im':float(S.imag)}

def run_case(path,levels=3200,seed=1729):
 t0=time.time();r=json.loads(path.read_text());P,Q=float(r['P']),float(r['Q'])
 w=float(r['positive_omega']);h=float(r['positive_eta']);Hs=float(r['S']['re']);target=wrap(float(r['phase_mod_2pi']))
 eng=PrefixProduct(levels);base=np.array([w,h-1j*np.pi,0],complex)
 def full(q,engine=eng):
  z=base+A@np.asarray(q);F,g,H=engine(z);S=F-z[0]*P-(z[1]+1j*np.pi)*Q
  return S,A.T@(g-np.array([P,Q,0])),A.T@H@A
 def gauge(u):return full([0,0,u])[1][2].real
 u=brentq(gauge,1e-5,.1,xtol=1e-15)
 q0=np.array([0,0,u]);S0,J0,K0=full(q0)
 eig,V=np.linalg.eigh(-K0.real)
 if np.min(eig)<=0:raise RuntimeError('positive peak Hessian is not negative definite')
 W=V@np.diag(1/np.sqrt(eig));cache={};calls=0
 def calc(p):
  nonlocal calls
  key=tuple(np.asarray(p))
  if key not in cache:
   S,J,K=full(q0+W@p);cache[key]=(S,W.T@J,W.T@K@W);calls+=1
  return cache[key]
 def loss(p):return float(S0.real-calc(p)[0].real)
 def gradloss(p):return -calc(p)[1].real
 def phase(p):return float((calc(p)[0]-S0).imag)
 def gradphase(p):return calc(p)[1].imag
 def gauge_positive(p):return float((q0+W@p)[2]-1e-8)
 gaugecon={'type':'ineq','fun':gauge_positive,'jac':lambda p:W[2]}
 rng=np.random.default_rng(seed+int(Q));directions=[]
 N=24
 for i in range(N):
  z=1-2*(i+.5)/N;az=np.pi*(3-np.sqrt(5))*i
  directions.append(np.array([np.sqrt(1-z*z)*np.cos(az),np.sqrt(1-z*z)*np.sin(az),z]))
 samples=[]
 for radius in [2.,5.,10.,20.,40.,70.]:
  for d in directions:
   p=radius*d
   if gauge_positive(p)>0:
    S,J,K=calc(p);samples.append({'p':p.tolist(),'phase':phase(p),'loss':loss(p),'height_margin':float(S.real-Hs)})
 # Find local unconstrained maxima from a small independent multistart set.
 peaks=[]
 for p in [np.zeros(3)]+[3*d for d in directions[::4]]:
  if gauge_positive(p)<=0:continue
  fit=minimize(loss,p,jac=gradloss,method='SLSQP',constraints=[gaugecon],bounds=[(-120,120)]*3,options={'maxiter':45,'ftol':1e-10})
  S,J,K=calc(fit.x);peaks.append({'physical_angles':(q0+W@fit.x).tolist(),'p':fit.x.tolist(),'S':encode(S),'success':bool(fit.success),'message':fit.message,'gradient_norm':float(np.linalg.norm(J.real))})
 bestpeak=max(peaks,key=lambda a:a['S']['re']);bestheight=bestpeak['S']['re']
 # The optimizer sees the continuous local phase, not a discontinuous wrapped equality.
 candidates=[];attempts=[]
 for branch in [target,target+2*np.pi,target-2*np.pi]:
  ranked=sorted(samples,key=lambda a:abs(a['phase']-branch)+.0003*max(a['loss'],0))
  # Different starts reduce accidental selection of a single connected branch.
  starts=[np.array(a['p']) for a in ranked[:3]]
  gp=gradphase(np.zeros(3));den=np.dot(gp,gp)
  if den>1e-10:
   linear=branch*gp/den
   if np.linalg.norm(linear)<120 and gauge_positive(linear)>0:starts.append(linear)
  for p in starts:
   con={'type':'eq','fun':lambda p,b=branch:phase(p)-b,'jac':gradphase}
   fit=minimize(loss,p,jac=gradloss,method='SLSQP',constraints=[gaugecon,con],bounds=[(-120,120)]*3,options={'maxiter':100,'ftol':3e-9})
   S,J,K=calc(fit.x);pr=phase(fit.x)-branch
   tangent=J.imag;mu=float(np.dot(J.real,tangent)/max(np.dot(tangent,tangent),1e-300))
   _,_,VT=np.linalg.svd(tangent.reshape(1,3));tang=VT[1:].T
   curv=np.linalg.eigvalsh(tang.T@(K.real-mu*K.imag)@tang)
   row={'physical_angles':(q0+W@fit.x).tolist(),'p':fit.x.tolist(),'S':encode(S),'phase_branch':float(branch),'phase_relative_peak':phase(fit.x),'phase_residual':float(pr),'phase_mod_2pi_residual':wrap(S.imag-float(r['phase_mod_2pi'])),'height_margin':float(S.real-Hs),'loss_from_best_peak':float(bestheight-S.real),'success':bool(fit.success),'message':fit.message,'iterations':int(fit.nit),'kkt_residual':float(np.linalg.norm(J.real-mu*J.imag)),'lagrange_multiplier':mu,'tangent_hessian_eigenvalues':curv.tolist(),'normal_coordinate_residual':[0.,0.,0.]}
   attempts.append(row)
   if abs(pr)<3e-6 and gauge_positive(fit.x)>0:candidates.append(row)
 winning=max(candidates,key=lambda a:a['S']['re']) if candidates else None
 # Repeat final candidates with a doubled product prefix. This is a convergence check.
 long=PrefixProduct(levels*2);checks=[]
 for label,c in [('positive_peak',bestpeak),('winning_phase_candidate',winning)]:
  if c is None:continue
  S,J,K=full(c['physical_angles'],long)
  checks.append({'label':label,'levels':levels*2,'S':encode(S),'action_change':float(abs(S-complex(c['S']['re'],c['S']['im']))),'phase_mod_2pi_residual':wrap(S.imag-float(r['phase_mod_2pi'])),'physical_angles':c['physical_angles']})
 return {'P':P,'Q':Q,'n':r['n'],'c':float(r['c']),'target_file':path.name,'saddle_height':Hs,'target_phase':target,'omega_radius_log':w,'eta_radius_log':h,'levels':levels,'center_physical_angles':q0.tolist(),'whitening':W.tolist(),'negative_height_hessian_eigenvalues':eig.tolist(),'positive_peak':bestpeak,'unconstrained_multistarts':peaks,'unconstrained_height_margin':float(bestheight-Hs),'phase_candidates':candidates,'phase_optimization_attempts':attempts,'winning_phase_candidate':winning,'samples':samples,'doubled_prefix_checks':checks,'engine_calls':calls,'seconds':time.time()-t0,'global_certification':False,'coverage':'24 sphere directions at radii 2,5,10,20,40,70 in Hessian-whitened coordinates; 7 local peak starts; at least 3 constrained starts on each of three local phase branches. Only the positive-u patch was searched; exact gauge inversion supplies negative u.'}

def main():
 p=argparse.ArgumentParser();p.add_argument('--targets',type=Path,default=OUT.parent/'targets');p.add_argument('--levels',type=int,default=3200);p.add_argument('--q',type=int,nargs='*');args=p.parse_args()
 paths=sorted(a for a in args.targets.glob('target_Q*.json') if re.fullmatch(r'target_Q[0-9]+',a.stem))
 if args.q:paths=[a for a in paths if int(a.stem.split('Q')[-1]) in args.q]
 if not paths:raise SystemExit('No target_Q*.json inputs found')
 rows=[]
 for a in paths:
  print('search',a.name,flush=True);row=run_case(a,args.levels);rows.append(row)
  (OUT/f'result_Q{int(row["Q"])}.json').write_text(json.dumps(row,indent=2))
  print('done',row['Q'],'margin',row['unconstrained_height_margin'],'phase margin',None if row['winning_phase_candidate'] is None else row['winning_phase_candidate']['height_margin'],'calls',row['engine_calls'],'seconds',row['seconds'],flush=True)
 (OUT/'results.json').write_text(json.dumps({'cases':rows,'global_certification':False},indent=2))
if __name__=='__main__':main()
