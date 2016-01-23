import math
import scipy.optimize as opt
M=[1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.50]
C=[0.383,0.382,0.375,0.366,0.356,0.346,0.337,0.328,0.32,0.313,0.288]
RHO=1.206
G=9.8
deltat=0.001
E=200000
def get_tau(y):
	if y<=9300:
		return 288.9-6.328*10**-3*y
	elif y<=12000:
		return 230.2-6.328*10**-3*(y-9300)+1.172*10**-6*(y-9300)**2
	else:
		return 221.5
def get_cs(tau):
	return (1.404*287*tau)**0.5
def get_ma(v,cs):
	return v/cs
def get_cxon(ma):
	if ma<M[0]:
		return C[0]
	if ma>M[-1]:
		return C[-1]
	for index,m in enumerate(M):
		if ma<m:
			left=ma-M[index-1]
			right=m-ma
			if left<right:
				return C[index-1]
			else:
				return C[index]
def get_i(H):
	return 2.9-1.373*H+0.32*H**2-0.0267*H**3
def get_v(E,m):
	return (2*E/m)**0.5
def get_cx(i,cxon):
	return i*cxon
def get_S(d):
	return math.pi*d**2/4
def get_t(y,v,theta,i,m,S,SD):
	sum_s=0
	count=0
	while sum_s<=SD:
		tau=get_tau(y)
		cs=get_cs(tau)
		ma=get_ma(v,cs)
		cxon=get_cxon(ma)
		cx=get_cx(i,cxon)
		v-=(RHO*v**2/(2*m)*S*cx-G*math.sin(theta))*deltat
		theta-=G*math.cos(theta)/v*deltat
		s=v*deltat
		sum_s+=s
		y+=s*math.sin(theta)
		count+=1
	return deltat*count
def cal_c(x):
	m=x[0]
	lam_n=x[1]
	lam_b=x[2]
	lam_c=x[3]
	H=lam_n+lam_b-0.3
	i=get_i(H)
	v=get_v(E,m)
	theta=65/360*2*3.14
	S=get_S(0.035)
	return get_t(0,v,theta,i,m,S,2000)
if __name__=='__main__':
	pass
	
