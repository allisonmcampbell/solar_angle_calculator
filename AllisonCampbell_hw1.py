from numpy import *
from operator import itemgetter

def getw(phi,long,tr,td,n,sign,pst):	
	w = []
	a = []
	g = []
	s = []
	for i in range(len(n)):
		ww = []
		aa = []
		gg = []
		ss = []
		b = tr*(n[i]-1.)*360/365
		E = 229.2*(7.5e-5 + 1.868e-3*cos(b) - 3.2077e-2*sin(b) - 1.4615e-2*cos(2*b) - 4.089e-2*sin(2*b))*15/60
		print 'E = ',E*60/15
		dec = 23.45*sin(tr*(284+n[i])*360/365)
		print 'dec = ',dec
		w_ew = td*arccos(tan(dec*tr)/tan(phi*tr))
		for j in range(len(pst)):
			ee = pst[j]+120 + E - long  
			print 'solar time = ',ee  
			ww.append(ee)
			ae,ge = getg(ee,w_ew,phi,dec,tr,td)
			aa.append(ae)
			gg.append(ge)
			#print b,E,dec,w_ew,ee,ae,ge
			shadow = getpoints(ge,ae,sign)
			ss.append(shadow)
		w.append(ww)
		a.append(aa)
		g.append(gg)
		s.append(ss)
	return(w,a,g,s)


def getg(w,w_ew,phi,dec,tr,td):
	if abs(w) <= w_ew:
		c1 =1.
	if abs(w) > w_ew:
		c1 = -1.
	if (phi - dec) >= 0:
		c2 = 1
	if (phi - dec) < 0:
		c2 = -1
	if (w >= 0):
		c3 = 1
	if (w < 0):
		c3 = -1
	theta_z = (arccos(cos(w*tr)*cos(dec*tr)*cos(phi*tr) + sin(dec*tr)*sin(phi*tr)))*td 
	g_s0 = (arcsin(sin(w*tr)*cos(dec*tr)/sin(theta_z*tr)))*td  
	g = c1*c2*g_s0 + c3*(1.-c1*c2)*90.
	a = 90.-theta_z	
	print 'theta z = ',theta_z
	return(a,g)
	
def getpoints(gamma,alpha,sign):
	tr=pi/180
	td=180/pi
	
	h0 = sign[0][1]/cos(tr*gamma)
	s00 = sign[0][0] - sign[0][1]*tan(tr*gamma)
	s01 = sign[0][1] - sign[0][1]
	s02 = sign[0][2] - h0*tan(tr*alpha)
	s0 = (s00,s01,s02)

	h1 = sign[1][1]/cos(tr*gamma)
	s10 = sign[1][0] - sign[1][1]*tan(tr*gamma)
	s11 = sign[1][1] - sign[1][1]
	s12 = sign[1][2] - h1*tan(tr*alpha)
	s1 = (s10,s11,s12)

	h2 = sign[2][1]/cos(tr*gamma)
	s20 = sign[2][0] - sign[2][1]*tan(tr*gamma)
	s21 = sign[2][1] - sign[2][1]
	s22 = sign[2][2] - h2*tan(tr*alpha)
	s2 = (s20,s21,s22)

	h3 = sign[3][1]/cos(tr*gamma)
	s30 = sign[3][0] - sign[3][1]*tan(tr*gamma)
	s31 = sign[3][1] - sign[3][1]
	s32 = sign[3][2] - h3*tan(tr*alpha)
	s3 = (s30,s31,s32)

	
	shadow = (s0,s1,s2,s3)
	return(shadow)

def calcarea(i,j,shadow,detector,n,pst):
	area = 0
	p0 = shadow[i][j][0]
	p1 = shadow[i][j][1]
	p2 = shadow[i][j][2]
	p3 = shadow[i][j][3]


	d0 = detector[0]
	d1 = detector[1]
	d2 = detector[2]
	d3 = detector[3]
	
	xmax=max(d0[0],d1[0],d2[0],d3[0])
	xmin=min(d0[0],d1[0],d2[0],d3[0])
	zmax=max(d0[2],d1[2],d2[2],d3[2])
	zmin=min(d0[2],d1[2],d2[2],d3[2])

	base = abs(p1[0] - p3[0])
	height = abs(p2[2] - p1[2])
	fullarea = .5*base*height	

	i0 = 0
	i1 = 0
	i2 = 0 
	i3 = 0
	
	if p0[0] >= xmin and p0[0] <= xmax: 
		if p0[2] >= zmin and p0[2] < zmax:
			i0 = 1

	if p1[0] >= xmin and p1[0] <= xmax: 
		if p1[2] >= zmin and p1[2] < zmax:
			i1 = 1

	if p2[0] >= xmin and p2[0] <= xmax: 
		if p2[2] >= zmin and p2[2] < zmax:
			i2 = 1

	if p3[0] >= xmin and p3[0] <= xmax: 
		if p3[2] >= zmin and p3[2] < zmax:
			i3 = 1

	if i0 == 0 and i1 == 0 and i2 == 0 and i3 == 0:
		area = 0
	
	if i0 == 0 and i1 == 1 and i2 == 0 and i3 == 0:
		#print n[i],pst[j],' p1 in',p2[2],p1[2],(p2[2]-p1[2])
		if (p2[2] > p1[2]):
			print 'hiya'
			if p3[2] > zmin:
				print 'ruh-roh raggy!'
			if p3[2] < zmin:		
				area = 0.5*(p1[0]-xmin)*(p1[0]-xmin)*(height/base) + (p1[0]-xmin)*(p1[2]-zmin)
				print n[i],pst[j], area, area/5.
		if (p2[2] == p1[2]):
			ztop = p1[2]
			zbottom = zmin
			dz =ztop-zbottom

			xleft = p1[0]
			xright = xmin
			dx = xleft-xright

			area = dx*dz
			print n[i],pst[j], area, area/5.



	if i0 == 0 and i1 == 0 and i2 == 1 and i3 == 0:
		#print n[i],pst[j],' p2 in'
		if (p2[2]-zmin) < height:
			area = 0.5*(p2[2]-zmin)*(p2[2]-zmin)*(base/height)
			print n[i],pst[j], area, area/5.

		if (p2[2]-zmin) >= height:
			area = fullarea - 0.5*(p0[0]-xmax)*(p0[0]-xmax)*(height/base) - (zmin-p2[2])*(xmax-p2[0])
			print n[i],pst[j], area, area/5.

	if i0 == 0 and i1 == 1 and i2 == 1 and i3 == 0:
		print n[i],pst[j],' p1 and p2 in'

		area = fullarea + (p1[2]-zmin)*(p1[0]-p2[0])

	if i0 == 0 and i1 == 1 and i2 == 1 and i3 == 1:
		print n[i],pst[j],' p1, p2, and p3 in'

	return(area)


	

########################################################################
phi = 40.87
long = 124.1
n = arange(259.,260.,1.)
#pst = arange(-30.,0,.5)
#n =[2]
pst=[.2505]
tr=pi/180. #to radians
td=180./pi #to degrees
# 1st orientation
#sign = ((-7.5,8,6.5),(-7.5,8,8.),(-10.5,8.,8.),(-10.5,8.,6.5))
# 2nd orientation
sign = ((-7.5,11.,6.5),(-7.5,11.,8.),(-7.5,8.,8.),(-7.5,8.,6.5))
detector =((0.,0.,4.), (0.,0.,6.), (-2.5,0.,6.), (-2.5,0.,4.))

	

w,a,g,s=getw(phi,long,tr,td,n,sign,pst)

for i in range(len(n)):
	for j in range(len(pst)):
		area=calcarea(i,j,s,detector,n,pst)
		

		





