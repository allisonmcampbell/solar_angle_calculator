from matplotlib.pyplot import *
from numpy import *
from operator import itemgetter

def getsunsetHA(dec,phi,tr,td):
	w_s = []
	for i in range(len(dec)):
		HA = td*arccos(-1*tan(dec[i]*tr)*tan(phi*tr))
		w_s.append(HA)
	return(w_s)

def getclearness(Hbar,dec,phi,w_s,n,G,tr,td):
	K_t = []
	print 'now printing ws, Hbar MJ, H0 MJ, Kt'
	for i in range(len(n)):
		##H [J/m^2/day]
		H = Hbar[i]*(3.6e6)
		##H0 units [J/m^2/day]
		H0 = (24*3600*G/pi)*(1+0.033*cos(tr*n[i]*360/365))*(cos(tr*phi)*cos(tr*dec[i])*sin(tr*w_s[i])+(pi*w_s[i]/180)*sin(tr*phi)*sin(tr*dec[i]))
		clearness = H/H0
		K_t.append(clearness)
		print H*1e-6
	return(K_t)

def getHdH(K_t,w_s,tr,td):
	Hd_H = []
	for i in range(len(w_s)):
		if w_s[i] <= 81.4:
			ratio = 1.391 - 3.56*K_t[i] + 4.189*K_t[i]**2 - 2.137*K_t[i]**3
		if w_s[i] > 81.4:
			ratio = 1.311 - 3.022*K_t[i] + 3.427*K_t[i]**2 - 1.821*K_t[i]**3
		Hd_H.append(ratio)
	return(Hd_H)

def getRbar(w_s,beta,gamma,phi,dec,Hd_H,Hbar,rho,tr,td):
	R = []
	H = []
	print 'now printing HdH, Rb, Ht'

	for i in range(len(w_s)):
		RR = []
		HH = []
		for j in range(len(gamma)):
			RRR = []
			HHH = []
			for k in range(len(beta)):
				A = cos(tr*beta[k]) + tan(tr*phi)*cos(tr*gamma[j])*sin(tr*beta[k])
				B = cos(tr*w_s[i])*cos(tr*beta[k]) + tan(tr*dec[i])*sin(tr*beta[k])*cos(tr*gamma[j])
				C = sin(tr*beta[k])*sin(tr*gamma[j])/cos(tr*phi)

				sr_test = td*arccos((A*B + C*(A**2 - B**2 + C**2)**0.5)/(A**2 + C**2))
				ss_test = td*arccos((A*B - C*(A**2 - B**2 + C**2)**0.5)/(A**2 + C**2))
				sr = min(w_s[i],sr_test)				
				ss = min(w_s[i],ss_test)
				w_sr = abs(sr)
				w_ss = -1*abs(ss)
				if (A > 0 and B > 0) or A>=B:
					w_sr = -1*w_sr
					w_ss = -1*w_ss
				
				a = 0.409 + 0.5016*sin(tr*(w_s[i]-60))
				b = 0.6609 - 0.4767*sin(tr*(w_s[i]-60))
				d = sin(tr*w_s[i]) - (pi*w_s[i]/180)*cos(tr*w_s[i])
				ap = a - Hd_H[i]

				##G1 = G(w_ss,w_sr)
				G1 = (1./(2.*d))*((b*A/2. - ap*B)*(w_ss-w_sr)*(pi/180.) + (ap*A-b*B)*(sin(tr*w_ss)-sin(tr*w_sr)) - ap*C*(cos(tr*w_ss)-cos(tr*w_sr)) + (b*A/2)*(sin(tr*w_ss)*cos(tr*w_ss)-sin(tr*w_sr)*cos(tr*w_sr)) + (b*C/2)*(sin(tr*w_ss)**2-sin(tr*w_sr)**2))
				##G2 = G(w_ss,-w_s)
				G2 = (1./(2.*d))*((b*A/2. - ap*B)*(w_ss+w_s[i])*(pi/180.) + (ap*A-b*B)*(sin(tr*w_ss)-sin(tr*-1*w_s[i])) - ap*C*(cos(tr*w_ss)-cos(tr*-1*w_s[i])) + (b*A/2)*(sin(tr*w_ss)*cos(tr*w_ss)-sin(tr*-1*w_s[i])*cos(tr*-1*w_s[i])) + (b*C/2)*(sin(tr*w_ss)**2-sin(tr*-1*w_s[i])**2))
				##G3 = G(w_s,w_sr)
				G3 = (1./(2.*d))*((b*A/2. - ap*B)*(w_s[i]-w_sr)*(pi/180.) + (ap*A-b*B)*(sin(tr*w_s[i])-sin(tr*w_sr)) - ap*C*(cos(tr*w_s[i])-cos(tr*w_sr)) + (b*A/2)*(sin(tr*w_s[i])*cos(tr*w_s[i])-sin(tr*w_sr)*cos(tr*w_sr)) + (b*C/2)*(sin(tr*w_s[i])**2-sin(tr*w_sr)**2))
				GG = G2+G3

				if w_ss >= w_sr:
					D = max(0,G1)	
				if w_sr > w_ss:
					print 'hiya!'
					D = max(0,GG)				

				R_bar = D + Hd_H[i]*((1+cos(tr*beta[k]))/2) + rho*((1-cos(tr*beta[k]))/2)
				#to MJ
				YO = Hbar[i]*(3.6)
				H_t = YO*R_bar
				RRR.append(R_bar)
				HHH.append(H_t)
			RR.append(RRR)
			HH.append(HHH)
		R.append(RR)
		H.append(HH)
	return(R,H)  
	#all of month, beta, gamma	

def plotgamma(Hbar_t,beta,gamma,n):
	print 'H = ',Hbar_t[0][0][0]
	HH = []
	for i in range(len(gamma)):
		H = 0
		for j in range(len(n)):
			H = H + Hbar_t[j][i][0]
		HH.append(H/12)	
		if i > 0:	
			if HH[i] > HH[i-1]:
				maxgamma = gamma[i]
				maxHH = HH[i]

	print 'max annual ave = ',max(HH),maxHH,maxgamma
		
	plot(gamma,HH,marker = '.',markerfacecolor='blue',markersize=1)
	xlabel('Azimuth Angle [degrees]')
	ylabel('Annual Ave Solar Energy [MJ/m^2]')
	text(55, 22, 'Max Annual Ave = 22.5', color='k')
	text(60, 21.7, 'Max Azimuth = 0', color='k')
	text(60, 21.4, 'Tilt Angle = 28.1', color='k')
	title('Effect of Azimuth Angle on Solar Collector, El Paso',fontsize=18, color='k')
	grid(True)
	show()


def plotbeta(Hbar_t,beta,gamma,n):
	HH = []
	for i in range(len(beta)):
		H = 0
		for j in range(len(n)):
			H = H + Hbar_t[j][0][i]
		HH.append(H/12)
		if i > 0:	
			if HH[i] > HH[i-1]:
				maxbeta = beta[i]
				maxHH = HH[i]

	print 'max annual ave = ',max(HH),maxHH,maxbeta


	plot(beta,HH,marker = '.',markerfacecolor='blue',markersize=1)
	xlabel('Tilt Angle [degrees]')
	ylabel('Annual Ave Solar Energy [MJ/m^2]')
	text(60, 22, 'Max Annual Ave = 22.5', color='k')
	text(60, 21.7, 'Max Tilt = 28.1', color='k')
	title('Effect of Tilt Angle on Solar Collector, Arcata',fontsize=18, color='k')
	grid(True)
	show()

def plotwinter(Hbar_t,beta,n):
	HH = []
	for i in range(len(beta)):			
		H = (Hbar_t[0][0][i] + Hbar_t[1][0][i] + Hbar_t[2][0][i] + Hbar_t[9][0][i] + Hbar_t[10][0][i] + Hbar_t[11][0][i])/6
		HH.append(H)
		print 'H = ',HH[i],i
		if i > 0:
			if HH[i] > HH[i-1]:
				maxbeta = beta[i]
				maxHH = HH[i]

	print 'max annual ave = ',max(HH),maxHH,maxbeta


	plot(beta,HH,marker = '.',markerfacecolor='blue',markersize=1)
	xlabel('Tilt Angle [degrees]')
	ylabel('Half-Annual Ave Solar Energy (MJ/m^2)')
	text(7, 12.6, 'Max Winter Ave = 12.6', color='k')
	text(7, 12.3, 'Max Tilt = 53.1', color='k')
	title('Effect of Tilt Angle, Winter Months, Arcata',fontsize=18, color='k')
	grid(True)
	show()

def plotsummer(Hbar_t,beta,n):
	HH = []
	for i in range(len(beta)):			
		H = (Hbar_t[3][0][i] + Hbar_t[4][0][i] + Hbar_t[5][0][i] + Hbar_t[6][0][i] + Hbar_t[7][0][i] + Hbar_t[7][0][i])/6
		HH.append(H)
		if i > 0:
			if HH[i] > HH[i-1]:
				maxbeta = beta[i]
				maxHH = HH[i]
			

	print 'max annual ave = ',max(HH),maxHH,maxbeta


	plot(beta,HH,marker = '.',markerfacecolor='blue',markersize=1)
	xlabel('Tilt Angle [degrees]')
	ylabel('Half-Annual Ave Solar Energy (MJ/m^2)')
	text(50, 25, 'Max Summer Ave = 26.3', color='k')
	text(50, 24, 'Max Tilt = 6.4', color='k')
	title('Effect of Tilt Angle, Summer Months, El Paso',fontsize=18, color='k')
	grid(True)
	show()
	

########################################################################
##arcata
#phi = 40.9
##el paso
phi = 31.8
###solar constant, [W/m^2]
G = 1367
n = (17,47,75,105,135,162,198,228,258,288,318,344)
dec = (-20.9,-13.0,-2.4,9.4,18.8,23.1,21.2,13.5,2.2,-9.6,-18.9,-23.0)
#gamma = arange(-150.,150.,1)
gamma = [0]
tilt = arange(0,90,.1)
#tilt = [19,20,1.]
tr=pi/180. #to radians
td=180./pi #to degrees
###ARCATA, [kWh/m^2/day]###
#Hbar = (1.8,2.5,3.6,5,5.8,6,5.9,5,4.4,3.1,2,1.6)
##el paso, tx 
Hbar = (3.5, 4.5, 5.9, 7.1, 7.8, 8.0, 7.4, 6.8, 5.9, 4.9, 3.8, 3.2, 5.7)
rho = 0.2


w_s = getsunsetHA(dec,phi,tr,td)

K_t = getclearness(Hbar,dec,phi,w_s,n,G,tr,td)

Hd_H = getHdH(K_t,w_s,tr,td)

Rbar,Hbar_t = getRbar(w_s,tilt,gamma,phi,dec,Hd_H,Hbar,rho,tr,td)

#plotgamma(Hbar_t,tilt,gamma,n)

#plotbeta(Hbar_t,tilt,gamma,n)

#plotwinter(Hbar_t,tilt,n)
#plotsummer(Hbar_t,tilt,n)

## [month][gamma][beta]
		





