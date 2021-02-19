#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Menu "Polaris"
	"pol_ARPES_sim",  Polaris_start()
End

Function radialfunct(n,l,Z)
	variable n,l,Z //Z原子番号
	string rad_s
	rad_s="R"+num2str(n)+num2str(l)
	
	make /N=200/O $rad_s
	wave Rnl=$rad_s
	setscale/P x,0,0.1,Rnl
	
	variable a0=1 //4*pi*ep0*hbar_Planck^2/(q_el^2*m_el) //(Borh 半径で規格化?)
	variable Anl=( (2*Z/(n*a0))^3 * gamma(n-l)/(2*n*gamma(n+l+1)) )^0.5
	//variable Anl=(2/(n*a0))^1.5
	 Rnl=Anl*( 2*Z*x/(n*a0) )^l*exp( -Z*x/(n*a0) )*laguerreA(n-l-1, 2*l+1, 2*Z*x/(n*a0))
	 //Rnl=( 2*r/(n*a0) )l*exp( -r/(n*a0) )*laguerreA(n-l-1, 2*l+1, 2*r/(n*a0))
	 display $rad_s
	 
end



Function wavefunct(n,l,m,H)
	variable n,l,m,H //H原子番号
	
	string l_s, m_s, wav,sph
	if(l==0)
		l_s="s"
		m_s=""
	elseif(l==1)
		l_s="p"
		if(m==1)
			m_s="x"
		elseif(m==0)
			m_s="z"
		elseif(m==-1)
			m_s="y"
		endif
	elseif(l==2)
		l_s="d"
		if(m==2)
			m_s="x2my2"
		elseif(m==1)
			m_s="xz"
		elseif(m==0)
			m_s="3z2mr2"
		elseif(m==-1)
			m_s="yz"
		elseif(m==-2)
			m_s="xy"
		endif
	elseif(l==3)
		l_s="f"
		if(m==3)
			m_s="3x2ymy3"	
		elseif(m==2)
			m_s="xyz"
		elseif(m==1)
			m_s="5z2m3r2"
		elseif(m==0)
			m_s="5z3m3r2z"
		elseif(m==-1)
			m_s="5z2xm3r2x"
		elseif(m==-2)
			m_s="x2zmy2z"
		elseif(m==-3)
			m_s="x3m3xy2"
		endif
	elseif(l==4)
		l_s="g"
	elseif(l==5)
		l_s="h"
	else
		l_s=num2str(l)
		m_s=num2str(m)
	endif
	
	wav="phi"+num2str(n)+l_s+m_s
	sph=num2str(n)+l_s+m_s
	
//	wave w=$wav
	Make/N=(200,200,200)/D/O $wav
	wave w=$wav
	setscale/p x,-20,0.2,w
	setscale/p y,-20,0.2,w
	setscale/p z,-20,0.2,w
	
	
	Make/N=(200,200,200)/D/O $sph
	wave realsp=$sph
	setscale/p x,-1.5,0.015,realsp
	setscale/p y,-1.5,0.015,realsp
	setscale/p z,-1.5,0.015,realsp

	variable a=1 // 4*pi*ep0*hbar_Planck^2/(q_el^2*m_el) //(Borh 半径で規格化?)
	variable Anl=(  (2*H/(n*a))^3 * gamma(n-l)/(2*n*gamma(n+l+1)) )^0.5
	//variable Anl=(2/(n*a0))^1.5
	 //w=-Anl*( 2*H*(r^2+y^2+z^2)^0.5/(n*a) )^l*exp( -H*(x^2+y^2+z^2)^0.5/(n*a) )*laguerreA(n-l-1, 2*l+1, 2*H*(x^2+y^2+z^2)^0.5/(n*a))//* real(sphericalHarmonics(2, 0, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))))//*0.25*(15/pi)^0.5*(x^2-y^2)/(x^2+y^2+z^2)
	 //Rnl=( 2*r/(n*a0) )l*exp( -r/(n*a0) )*laguerreA(n-l-1, 2*l+1, 2*r/(n*a0))
	 
	 if(m<0)//好きな関数を作る絶対値など、実数球面調和関数は等直面プロットを使う。
	 	w=-Anl*( 2*H*(r^2+y^2+z^2)^0.5/(n*a) )^l*exp( -H*(x^2+y^2+z^2)^0.5/(n*a) )*laguerreA(n-l-1, 2*l+1, 2*H*(x^2+y^2+z^2)^0.5/(n*a)) *( -imag( (-1)^(-m) *sphericalHarmonics(l, abs(m), acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))-sphericalHarmonics(l, -abs(m), acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5)) ) ) 
	 	realsp= -imag( (-1)^(-m) *sphericalHarmonics(l, abs(m), acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))-sphericalHarmonics(l, -abs(m), acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5)) ) /(x^2+y^2+z^2)^0.5
	 elseif(m==0)
	 	w=-Anl*( 2*H*(r^2+y^2+z^2)^0.5/(n*a) )^l*exp( -H*(x^2+y^2+z^2)^0.5/(n*a) )*laguerreA(n-l-1, 2*l+1, 2*H*(x^2+y^2+z^2)^0.5/(n*a)) *real(sphericalHarmonics(l, m, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5)))
	 	realsp=real(sphericalHarmonics(l, m, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))) /(x^2+y^2+z^2)^0.5
	 elseif(m>0)
	 	w=-Anl*( 2*H*(r^2+y^2+z^2)^0.5/(n*a) )^l*exp( -H*(x^2+y^2+z^2)^0.5/(n*a) )*laguerreA(n-l-1, 2*l+1, 2*H*(x^2+y^2+z^2)^0.5/(n*a)) *( real( (-1)^(-m) *sphericalHarmonics(l, abs(m), acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))+sphericalHarmonics(l, -abs(m), acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5)) ) ) 
	 	realsp=real( (-1)^(-m) *sphericalHarmonics(l, abs(m), acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))+sphericalHarmonics(l, -abs(m), acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5)) )/(x^2+y^2+z^2)^0.5
	 endif
end

function Gegenbaur(m,v,wav)
	variable m,v
	string wav
		Make/N=(200,200)/D/O $wav
	wave w=$wav
	setscale/p x,-20,0.2,w
	setscale/p y,-20,0.2,w
	//variable G
		w=gamma(m+1/2)*gamma(v+2*m)/(gamma(v+1)*gamma(2*m))* ((1-y^2)/4)^(0.25-0.5*m)*legendreA(0.5-m, v+m-0.5, y )
		

end


function crosssection(n,l,wav)
	variable n,l
	string wav
	Make/N=1000/D/O $wav
	wave w=$wav
	setscale/p x,4.8,0.02,w
	variable a=0.529177210903
	
	w=(n*a*0.5123*(x-4.8)^0.5)^(2*l)/((n*a*0.5123*(x-4.8)^0.5)^2+1)^(2*l+4) * (gamma(l+1+1/2)*gamma(n-l-1+2*(l+1))/(gamma(n-l)*gamma(2*l+2)))^2   * (   ((1-   (((n*a*0.5123*(x-4.8)^0.5)^2-1)/((n*a*0.5123*(x-4.8)^0.5)^2+1))^2)/4)^(0.25-0.5*(l+1))*legendreA(n-l-1+(l+1)-0.5, 0.5-l-1,((n*a*0.5123*(x-4.8)^0.5)^2-1)/((n*a*0.5123*(x-4.8)^0.5)^2+1) ))^2


end


///                 (  ( (n*a*y)^2-1 )/( (n*a*y)^2+1) )^2 )

Function projection_sph_kx_ky_DA30(l,m,hn,W) //angular oribital term
	variable l,m,hn,W 
	variable kf=0.512*(hn-W)^0.5
	string l_s, m_s, sph
	if(l==0)
		l_s="s"
		m_s=""
	elseif(l==1)
		l_s="p"
		if(m==1)
			m_s="x"
		elseif(m==0)
			m_s="z"
		elseif(m==-1)
			m_s="y"
		endif
	elseif(l==2)
		l_s="d"
		if(m==2)
			m_s="x2my2"
		elseif(m==1)
			m_s="xz"
		elseif(m==0)
			m_s="3z2mr2"
		elseif(m==-1)
			m_s="yz"
		elseif(m==-2)
			m_s="xy"
		endif
	elseif(l==3)
		l_s="f"
		if(m==3)
			m_s="3x2ymy3"	
		elseif(m==2)
			m_s="xyz"
		elseif(m==1)
			m_s="5z2m3r2"
		elseif(m==0)
			m_s="5z3m3r2z"
		elseif(m==-1)
			m_s="5z2xm3r2x"
		elseif(m==-2)
			m_s="x2zmy2z"
		elseif(m==-3)
			m_s="x3m3xy2"
		endif
	elseif(l==4)
		l_s="g"
	elseif(l==5)
		l_s="h"
	else
		l_s=num2str(l)
		m_s=num2str(m)
	endif
	
	sph=l_s+m_s
	
	Make/N=(200,200)/D/O $sph
	wave realsp=$sph
	setscale/p x,-1.5,0.015,realsp ////あとで編集
	setscale/p y,-1.5,0.015,realsp ///あとで編集

	//m=-2,-1,,0,1,2,立方調和関数で定義。例m=2の時dx^2-y^2
	
	 
	 if(m<0)//好きな関数を作る絶対値など、実数球面調和関数は等直面プロットを使う。
	 	 
	 	realsp= abs(-imag( (-1)^(-m) *sphericalHarmonics(l, abs(m), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5))-sphericalHarmonics(l, -abs(m), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5)) )/2^0.5) 
	 elseif(m==0)
	 	//(1-x^-y^2)^0.5
	 	//realsp=abs(real(sphericalHarmonics(l, m, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))))*realsp
	 	realsp=abs(real(sphericalHarmonics(l, m, asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5))))
	 elseif(m>0)
	 	
	 	realsp=abs(real( (-1)^(-m) *sphericalHarmonics(l, abs(m), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5))+sphericalHarmonics(l, -abs(m), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5)) )/2^0.5)
	 endif
	 
	 //for
	 //endfor
end

///////////偏光ARPESmacro開発中////////////////
////グローバル変数作成パート//////////////
Macro Polaris_start()
	if(strlen(WinList("POLARIS_window", ";", "WIN:64")))
			DoWindow/F POLARIS_window
	else
		String nf=GetDataFolder(1)
		NewDataFolder/O root:Polaris
		NewDataFolder/O root:Polaris:misc
		NewDataFolder/O root:Polaris:outuput
		NewDataFolder/O root:Polaris:misc:light_pol
		SetDataFolder root:Polaris:misc:light_pol
		variable/g v_guzal=0
		variable/g v_Delta=0
		variable/g v_inc_angle=50 //(deg)
		variable/g v_inc_tilt_angle=0
		
		
		//NewDataFolder/O root:pol_ARPES:misc:rot_matrix
		//NewDataFolder/O root:pol_ARPES:Curves
		NewDataFolder/S/O root:Polaris:global
		//make/T/O/N=0 $("win_list")
		SetDataFolder root:Polaris:misc
		variable/g v_hn
		if(v_hn==0)
			v_hn=50
		endif
		variable/g v_V0
		if(v_V0==0)
			v_V0=12
		endif
		variable/g v_W
		if(v_W==0)
			v_W=4.5
		endif
		
		string/g st_orb
		
		variable/g v_escdep=10
		variable/g v_escdep_fit=(143/(50-4.5)^2+0.054*(50-4.5)^0.5)*10
		variable/g v_Intofs=0
		//variable/g v_EB
		variable/g v_theta_s
		variable/g v_theta_e=-15+0.1*(301-1)
		if(v_theta_s==0)
			variable/g v_theta_s=-15
			//variable/g v_theta_e=15
		endif
		variable/g v_tilt_s
		variable/g v_tilt_e=-15+0.1*(301-1)
		if(v_tilt_s==0)
			variable/g v_tilt_s=-15
			//variable/g v_theta_e=15
		endif
		variable/g v_tilt=0
		
		//wave作成用変数
		variable/g v_kind=2
		
		variable/g v_Epoint=101
		variable/g v_Edel=0.05
		variable/g v_Est=-2
		variable/g v_Een=-2+(101-1)*0.05
		
		variable/g v_tilt_point=301
		variable/g v_tilt_del=0.1
		
		variable/g v_polar_point=301
		variable/g v_polar_del=0.1
		
		variable/g v_azi=0
		//angle or k mode切り替え
		variable/g v_check_angle_or_k=0
		
		//variable/g v_pol_ofs=0
		//variable/g v_tilt_ofs=0
		//variable/g v_azi_ofs=0
		
		variable/g v_l=1
		variable/g v_s=1
		
		//wave作成用変数
		variable/g v_kind=1
		
		variable/g v_Epoint=100
		variable/g v_Edel=0.05
		variable/g v_Est=-2
		
		variable/g v_kpoint=180
		variable/g v_kdel=1
		variable/g v_kst=-90
		
		variable/g tiltmap2DA30=0
		variable/g v_liveupdate=0
		variable/g v_select_display=3
		SetDataFolder root:Polaris:outuput
		//PrepareWave(v_kind, v_Epoint,v_Est,v_Edel, v_polar_point, v_theta_s, v_polar_del,v_tilt_point, v_tilt_s, v_tilt_del)
		//wave pre_w
		//duplicate/o pre_w, polwav
		//duplicate/o pre_w, $cubspfunc
		//duplicate/o pre_w, $tot
		POLARIS_MACRO()
		
	endif

End
	



////偏光パート/////////
function Polaris_pol_light_DA30(Lpol, Cpol, incdeg_pol,incdeg_til, V0,escdep, hn,W,x,y,z,kind, check_angle_or_k)
	variable Lpol,Cpol,incdeg_pol,incdeg_til,V0, escdep, hn,W,x,y,z,kind,check_angle_or_k ///偏光角度(deg)、入射角、光電子脱出深さ、内部ポテンシャル、励起光、仕事関数
	//Lpol=pi/2: p-pol, スリットと平行、Lpol=0:s-pol
	//Lpol=pi/4かつCpol=pi/2:right cir_pol.
	//Lpol=pi/4かつCpol=-pi/2:left cir_pol. 
	variable kf=0.512*(hn-W)^0.5
	variable polwav
	//使うのは p-pol, s-pol, cir_polm 逆
	//変数:入射光の角度、光電子の脱出深さ
	

	//NVAR check_angle_or_k=root:pol_ARPES:misc:v_check_angle_or_k
	 //variable check_angle_or_k=1
	 //variable kind=2
	//x=kfsin(pol)
	//y=kfcospol cos tilt
	//z=kfcospol sin tilt
	//////関数代入パート
	variable S0,S1,S2,S3, E_int
	E_int=1 //とりま1
	variable a, lamda
	a=incdeg_pol; lamda=escdep
	S0=E_int^2
	S1=E_int^2*cos(2*Lpol*pi/180)
	S2=E_int^2*cos(Cpol*pi/180)*sin(2*Lpol*pi/180) // 円偏光
	S3=E_int^2*sin(Cpol*pi/180)*sin(2*Lpol*pi/180) //円偏光
	
	variable t=incdeg_til
	variable kf_x=0.512*(hn-W+x)^0.5
	if(check_angle_or_k==0)///DA30angle
		if(kind==1)
		//kind1E_k
	 		variable theta_Ek=((y*pi/180)^2+(0*pi/180)^2)^0.5
	 		variable psi_Ek=Polaris_nan2zero(sign(sin(0*pi/180))*acos(y/(y^2+0^2)^0.5*sin((y^2+0^2)^0.5*pi/180)/(y^2/(y^2+0^2)*sin((y^2+0^2)^0.5*pi/180)^2+0^2/(y^2+0^2)*sin((y^2+0^2)^0.5*pi/180)^2)^0.5))
	 		variable kx_E=kf_x*sin(theta_Ek)*cos(psi_Ek)
	 		variable ky_E=kf_x*sin(theta_Ek)*sin(psi_Ek)
	 		variable kz_E=kf_x*cos(theta_Ek)
			polwav=S0/2*((0*cos(t*pi/180))^2+  (kx_E*cos(pi/180*a)+(kf_x^2-(kx_E^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((0*cos(t*pi/180))^2-  (kx_E*cos(pi/180*a)+(kf_x^2-(kx_E^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(0*cos(t*pi/180))*(kx_E*cos(pi/180*a)+(kf_x^2-(kx_E^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(0*cos(t*pi/180))*sin(pi/180*a)/lamda
		elseif(kind==2)
		//kind2 kx-ky
			variable theta_k=((x*pi/180)^2+(y*pi/180)^2)^0.5
	 		variable psi_k=Polaris_nan2zero(sign(sin(y*pi/180))*acos(x/(x^2+y^2)^0.5*sin((x^2+y^2)^0.5*pi/180)/(x^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2+y^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2)^0.5))
			variable kx=kf*sin(theta_k)*cos(psi_k)
	 		variable ky=kf*sin(theta_k)*sin(psi_k)
	 		variable kz=kf*cos(theta_k)
			polwav=S0/2*((ky*cos(t*pi/180))^2+  (kx*cos(pi/180*a)+(kf^2-(kx^2+ky^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((ky*cos(t*pi/180))^2-  (kx*cos(pi/180*a)+(kf^2-(kx^2+ky^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(ky*cos(t*pi/180))*(kx*cos(pi/180*a)+(kf^2-(kx^2+ky^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(ky*cos(t*pi/180))*sin(pi/180*a)/lamda
		elseif(kind==3)
		//E(kx,ky)
			variable theta_k_vol=((y*pi/180)^2+(z*pi/180)^2)^0.5
	 		variable psi_k_vol=Polaris_nan2zero(sign(sin(z*pi/180))*acos(y/(y^2+z^2)^0.5*sin((y^2+z^2)^0.5*pi/180)/(y^2/(y^2+z^2)*sin((y^2+z^2)^0.5*pi/180)^2+z^2/(y^2+z^2)*sin((y^2+z^2)^0.5*pi/180)^2)^0.5))
			variable kx_vol=kf*sin(theta_k_vol)*cos(psi_k_vol)
	 		variable ky_vol=kf*sin(theta_k_vol)*sin(psi_k_vol)
	 		variable kz_vol=kf*cos(theta_k_vol)
			polwav=S0/2*((ky_vol*cos(t*pi/180))^2+  (kx_vol*cos(pi/180*a)+(kf_x^2-(kx_vol^2+ky_vol^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((ky_vol*cos(t*pi/180))^2-  (kx_vol*cos(pi/180*a)+(kf_x^2-(kx_vol^2+ky_vol^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(ky_vol*cos(t*pi/180))*(kx_vol*cos(pi/180*a)+(kf_x^2-(kx_vol^2+ky_vol^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(ky_vol*cos(t*pi/180))*sin(pi/180*a)/lamda
		else
			abort
		endif
		
	elseif(check_angle_or_k==1)//DA30 k
		if(kind==1)
			polwav=S0/2*((0*cos(t*pi/180))^2+  (y*cos(pi/180*a)+(kf_x^2-(y^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((0*cos(t*pi/180))^2-  (y*cos(pi/180*a)+(kf_x^2-(y^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(0*cos(t*pi/180))*(y*cos(pi/180*a)+(kf_x^2-(y^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(0*cos(t*pi/180))*sin(pi/180*a)/lamda
		elseif(kind==2)
			polwav=S0/2*((y*cos(t*pi/180))^2+  (x*cos(pi/180*a)+(kf^2-(x^2+y^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((y*cos(t*pi/180))^2-  (x*cos(pi/180*a)+(kf^2-(x^2+y^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(y*cos(t*pi/180))*(x*cos(pi/180*a)+(kf^2-(x^2+y^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(y*cos(t*pi/180))*sin(pi/180*a)/lamda
		elseif(kind==3)
			//variable kf_vol=0.512*(hn-W+x)^0.5
			polwav=S0/2*((z*cos(t*pi/180))^2+  (y*cos(pi/180*a)+(kf_x^2-(y^2+z^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((z*cos(t*pi/180))^2-  (y*cos(pi/180*a)+(kf_x^2-(y^2+z^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(z*cos(t*pi/180))*(y*cos(pi/180*a)+(kf_x^2-(y^2+z^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(z*cos(t*pi/180))*sin(pi/180*a)/lamda
		else
			abort
		endif
	//polwav=-S3*x*sin(pi/180*a)/lamda+S0/2*(x^2+  (y*cos(pi/180*a)+0*(kf-(x^2+y^2)^0.5+kf*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2)
	endif
	return polwav
	//polwav=S0/2*(x^2+  (y*cos(pi/180*a)+(kf-(x^2+y^2)^0.5+kf*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) + S1/2*(x^2-  (y*cos(pi/180*a)+((kf-(x^2+y^2)^0.5+kf*V0)^0.5)*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) +S2*x*(y*cos(pi/180*a)+(kf-(x^2+y^2)^0.5+kf*V0)^0.5*sin(pi/180*a))-S3*x*sin(pi/180*a)/lamda
End



function Polaris_pol_light_tiltmap(Lpol, Cpol, incdeg_pol, V0,escdep, hn,W,x,y,z,kind,check_angle_or_k)
	variable Lpol,Cpol,incdeg_pol, V0,escdep, hn,W,x,y,z,kind,check_angle_or_k ///偏光角度(deg)、入射角、光電子脱出深さ、内部ポテンシャル、励起光、仕事関数
	//Lpol=pi/2: p-pol, スリットと平行、Lpol=0:s-pol
	//Lpol=pi/4かつCpol=pi/2:right cir_pol.
	//Lpol=pi/4かつCpol=-pi/2:left cir_pol. 
	variable kf=0.512*(hn-W)^0.5
	variable polwav

	//使うのは p-pol, s-pol, cir_polm 逆
	//変数:入射光の角度、光電子の脱出深さ
	

	//NVAR check_angle_or_k=root:pol_ARPES:misc:v_check_angle_or_k
	 //variable check_angle_or_k=1
	 //variable kind=2
	//x=kfsin(pol)
	//y=kfcospol cos tilt
	//z=kfcospol sin tilt
	//////関数代入パート
	variable S0,S1,S2,S3, E_int
	E_int=1 //とりま1
	variable a, lamda
	a=incdeg_pol; lamda=escdep
	S0=E_int^2
	S1=E_int^2*cos(2*Lpol*pi/180)
	S2=E_int^2*cos(Cpol*pi/180)*sin(2*Lpol*pi/180) // 円偏光
	S3=E_int^2*sin(Cpol*pi/180)*sin(2*Lpol*pi/180) //円偏光
	
	///variable t=incdeg_til
	variable kf_x=0.512*(hn-W+x)^0.5
	if(check_angle_or_k==0)///DA30angle
		if(kind==1)
		//kind1E_k
	 		variable theta_Ek=((y*pi/180)^2+(0*pi/180)^2)^0.5
	 		variable psi_Ek=Polaris_nan2zero(sign(sin(0*pi/180))*acos(y/(y^2+0^2)^0.5*sin((y^2+0^2)^0.5*pi/180)/(y^2/(y^2+0^2)*sin((y^2+0^2)^0.5*pi/180)^2+0^2/(y^2+0^2)*sin((y^2+0^2)^0.5*pi/180)^2)^0.5))
	 		variable kx_E=kf_x*sin(theta_Ek)*cos(psi_Ek)
	 		variable ky_E=kf_x*sin(theta_Ek)*sin(psi_Ek)
	 		variable kz_E=kf_x*cos(theta_Ek)
			polwav=S0/2*((0*cos(0*pi/180))^2+  (kx_E*cos(pi/180*a)+(kf_x^2-(kx_E^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((0*cos(0*pi/180))^2-  (kx_E*cos(pi/180*a)+(kf_x^2-(kx_E^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(0*cos(0*pi/180))*(kx_E*cos(pi/180*a)+(kf_x^2-(kx_E^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(0*cos(0*pi/180))*sin(pi/180*a)/lamda
		elseif(kind==2)//書き換え終了
		//kind2 kx-ky
			variable theta_k=((x*pi/180)^2+(-y*pi/180)^2)^0.5
	 		variable psi_k=Polaris_nan2zero(sign(sin(-y*pi/180))*acos(x/(x^2+y^2)^0.5*sin((x^2+y^2)^0.5*pi/180)/(x^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2+y^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2)^0.5))
			variable kx=kf*sin(theta_k)*cos(psi_k)
	 		variable ky=kf*sin(theta_k)*sin(psi_k)
	 		variable kz=kf*cos(theta_k)
			polwav=S0/2*((ky*cos(-y*pi/180))^2+  (kx*cos(pi/180*a)+(kf^2-(kx^2+ky^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((ky*cos(y*pi/180))^2-  (kx*cos(pi/180*a)+(kf^2-(kx^2+ky^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(ky*cos(-y*pi/180))*(kx*cos(pi/180*a)+(kf^2-(kx^2+ky^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(ky*cos(-y*pi/180))*sin(pi/180*a)/lamda
		elseif(kind==3)
		//E(kx,ky)
			variable theta_k_vol=((y*pi/180)^2+(-z*pi/180)^2)^0.5
	 		variable psi_k_vol=Polaris_nan2zero(sign(sin(-z*pi/180))*acos(y/(y^2+z^2)^0.5*sin((y^2+z^2)^0.5*pi/180)/(y^2/(y^2+z^2)*sin((y^2+z^2)^0.5*pi/180)^2+z^2/(y^2+z^2)*sin((y^2+z^2)^0.5*pi/180)^2)^0.5))
			variable kx_vol=kf*sin(theta_k_vol)*cos(psi_k_vol)
	 		variable ky_vol=kf*sin(theta_k_vol)*sin(psi_k_vol)
	 		variable kz_vol=kf*cos(theta_k_vol)
			polwav=S0/2*((ky_vol*cos(-z*pi/180))^2+  (kx_vol*cos(pi/180*a)+(kf_x^2-(kx_vol^2+ky_vol^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((ky_vol*cos(-z*pi/180))^2-  (kx_vol*cos(pi/180*a)+(kf_x^2-(kx_vol^2+ky_vol^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(ky_vol*cos(-z*pi/180))*(kx_vol*cos(pi/180*a)+(kf_x^2-(kx_vol^2+ky_vol^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(ky_vol*cos(-z*pi/180))*sin(pi/180*a)/lamda
		else
			abort
		endif
		
	elseif(check_angle_or_k==1)//DA30 k
		if(kind==1)//G-lineのみ
			polwav=S0/2*((0*cos(0*pi/180))^2+  (y*cos(pi/180*a)+(kf_x^2-(y^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((0*cos(0*pi/180))^2-  (y*cos(pi/180*a)+(kf_x^2-(y^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(0*cos(0*pi/180))*(y*cos(pi/180*a)+(kf_x^2-(y^2+0^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(0*cos(0*pi/180))*sin(pi/180*a)/lamda
		elseif(kind==2)
			variable t_kxky=180/pi*sign(-y)*asin(-y/(kf^2-x^2)^0.5)
			
			polwav=S0/2*((-y*cos(t_kxky*pi/180))^2+  (x*cos(pi/180*a)+(kf^2-(x^2+y^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((-y*cos(t_kxky*pi/180))^2-  (x*cos(pi/180*a)+(kf^2-(x^2+y^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(-y*cos(t_kxky*pi/180))*(x*cos(pi/180*a)+(kf^2-(x^2+y^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(-y*cos(t_kxky*pi/180))*sin(pi/180*a)/lamda
		elseif(kind==3)
			//variable kf_vol=0.512*(hn-W+x)^0.5
			variable t_vol=180/pi*sign(y)*asin(-z/(kf_x^2-y^2)^0.5)
			polwav=S0/2*((-z*cos(t_vol*pi/180))^2+  (y*cos(pi/180*a)+(kf_x^2-(y^2+z^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) \
			+S1/2*((-z*cos(t_vol*pi/180))^2-  (y*cos(pi/180*a)+(kf_x^2-(y^2+z^2)+0.512^2*V0)^0.5*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) \
			+S2*(-z*cos(t_vol*pi/180))*(y*cos(pi/180*a)+(kf_x^2-(y^2+z^2)+0.512^2*V0)^0.5*sin(pi/180*a))\
			-S3*(-z*cos(t_vol*pi/180))*sin(pi/180*a)/lamda
		else
			abort
		endif
	//polwav=-S3*x*sin(pi/180*a)/lamda+S0/2*(x^2+  (y*cos(pi/180*a)+0*(kf-(x^2+y^2)^0.5+kf*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2)
	endif
	return polwav
	//polwav=S0/2*(x^2+  (y*cos(pi/180*a)+(kf-(x^2+y^2)^0.5+kf*V0)^0.5*sin(pi/180*a))^2 +(sin(pi/180*a)/lamda)^2) + S1/2*(x^2-  (y*cos(pi/180*a)+((kf-(x^2+y^2)^0.5+kf*V0)^0.5)*sin(pi/180*a))^2 -(sin(pi/180*a)/lamda)^2) +S2*x*(y*cos(pi/180*a)+(kf-(x^2+y^2)^0.5+kf*V0)^0.5*sin(pi/180*a))-S3*x*sin(pi/180*a)/lamda
End




Function Polaris_Wigner_Rot_Matrix_DA30(v_l,v_m,v_hn,v_W, v_pol,v_tilt,v_azi,check_angle_or_k, kind, x,y,z)///実装用
	variable/D v_l,v_m,v_hn,v_W, v_pol,v_tilt,v_azi,check_angle_or_k,kind, x,y,z
	variable/D realsp
	v_tilt=(-1)*v_tilt
	
	variable/D kf=0.512*(v_hn-v_W)^0.5//(エネルギー方向の変数を入れる)
	variable/D a, b, c, d, e, f, g, h, o //回転行列の要素
	//a= mw_tot[0][0]; b= mw_tot[0][1]; c= mw_tot[0][2]
	//d= mw_tot[1][0]; e= mw_tot[1][1]; f= mw_tot[1][2]
	//g= mw_tot[2][0]; h= mw_tot[2][1]; o= mw_tot[2][2]
	
 	a= Cos(v_azi*pi/180)*Cos(v_pol*pi/180) - Sin(v_azi*pi/180)*Sin(v_pol*pi/180)*Sin((v_tilt)*pi/180)
	b= -Cos((v_tilt)*pi/180)* Sin(v_azi*pi/180)
	c= Cos(v_azi*pi/180)* Sin(v_pol*pi/180) + Cos(v_pol*pi/180)* Sin(v_azi*pi/180)* Sin((v_tilt)*pi/180)
	d= Cos(v_pol*pi/180)* Sin(v_azi*pi/180) + Cos(v_azi*pi/180)* Sin(v_pol*pi/180)* Sin((v_tilt)*pi/180)
	
	e= Cos(v_azi*pi/180)* Cos((v_tilt)*pi/180) 
	f= Sin(v_azi*pi/180)* Sin(v_pol*pi/180) - Cos(v_azi*pi/180)* Cos(v_pol*pi/180)* Sin((v_tilt)*pi/180)
	g= -Cos((v_tilt)*pi/180)* Sin(v_pol*pi/180)
	h= Sin((v_tilt)*pi/180)
	o= Cos(v_pol*pi/180)* Cos((v_tilt)*pi/180)
	
	//Ttot=R(azi)R(tilt)R(polar) (回転行列の順番はマニピュレータ依存。HiSORの場合tiltとaziはおそらく逆)sampleの傾きの補正も必要なら6つの回転行列が必要になる。
	//Ttot={(a,b,c,), (d,e,f),(g,h,0)
	//将来的に回転要素をxyz好きに組めると良いかも
	//Wigner D Matrix
	//Rem. only s-p-d- orbital. hard to express f-orbital. if you want to use f-orbital,please write down.
	//Rem. about (l,m), 
	//Def.	  (0,0)=s, 
			//(1,1)=(px), (1,0)=(pz), (1,-1)=(py)
			//(2,2)=(d_x^2-y^2), (2,1)=(d_zx),(2,0)=(d_z^2), (2,-1)=(d_yz), (2,-2)=(d_xy)
			

	
	//NVAR check_angle_or_k=root:pol_ARPES:misc:v_check_angle_or_k
	

	/////////////関数代入パート///////////////////
	if(v_l==0 && v_m==0) //s-orb nochange
		if(kind==1)
			realsp=abs(PrepareRealSpH(0,0,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k) )
		elseif(kind==2)
			realsp=abs(PrepareRealSpH(0,0,kf,x,y,check_angle_or_k) )
		elseif(kind==3)
			realsp=abs(PrepareRealSpH(0,0,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k) )
		else
			Abort
		endif
	elseif(v_l==1 && -abs(v_l)<=v_m && v_m<=v_l) ///p-orb part
     //////px,py, pz////選択パート
		variable/D t, u, v
		if(v_m==1)//px
			t=a; u=d; v=g
		elseif(v_m==-1)//py
			t=b; u=e; v=h
		else //pz
			t=c; u=f; v=o
		endif
	///////////////////sample Rot. part////////////
	///tilymap時サンプルに対する検出面が変わってしまう問題をどうするか? =>for分で回す(脳筋)
	//		解決方法1. tiltabgele をfor文を回して連続変化させ、その都度ky(y)=kf*cos(pi/180*polar)*sin(pi/180*tilt)を代入:polar=constとする
	//ただし、データ点の間隔が一定で無くなる。=>kyを連続変化させれば良い設定すればよくない？
	
	//解決方法2. 計算して解析的に求める///むずい
		//realsp=abs(1/3*(t*(real( (-1)^(-1) *sphericalHarmonics(l, abs(1), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5))+sphericalHarmonics(l, -abs(1), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5)) )/2^0.5) \
				//+v*(real(sphericalHarmonics(l, 0, asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5)))) \
				//+u*(-imag( (-1)^(-(-1)) *sphericalHarmonics(l, abs(-1), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5))-sphericalHarmonics(l, -abs(-1), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5)) )/2^0.5) )\
				//)
				
	////////////if文にてkx-kyかE-kかvolかtiltmapかを判断			
			
	switch(kind)
		case 1:	//E-kx
			realsp=abs((t*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+v*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+u*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)))
		break
			
		case 2://only FS
		
			realsp=abs((t*PrepareRealSpH(v_l,1,kf,x,y,check_angle_or_k)+v*PrepareRealSpH(v_l,0,kf,x,y,check_angle_or_k)+u*PrepareRealSpH(v_l,-1,kf,x,y,check_angle_or_k)))
			//realsp=abs(PrepareRealSpH(v_l,0,kf,x,y,check_angle_or_k))
			
			//realsp=abs(1/3*(t*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+v*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+u*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)))
		break
		case 3://3Dvol
			realsp=abs(1/3*(t*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W-x)^0.5,y,z,check_angle_or_k)+v*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W-x)^0.5,y,z,check_angle_or_k)+u*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W-x)^0.5,y,z,check_angle_or_k)))
		break
	endswitch
		
		
////////d-orbパート//////////////////////////	
	elseif(v_l==2 && -abs(v_l)<=v_m && v_m<=v_l) ///d-orb part
		variable/D d1, d2, d3, d4, d5
		if(v_m==2)//dx^2-y^2 rot
			d1=a*d-b*e //dxy'
			d2=0.5*(a^2-b^2-d^2+e^2) //dx^2-y^2'
			d3=a*g-b*h //dzx'
			d4=d*g-e*h //dyz'
			d5=-3^0.5/6*(a^2-b^2+d^2-e^2-2*g^2+2*h^2) //dz^2'
			
		elseif(v_m==1)//dzx
			d1=c*d+a*f //dxy'
			d2=a*c-d*f
			d3=c*g+a*o
			d4=f*g+d*o
			d5=-3^0.5/3*(a*c+d*f-2*g*o)
		elseif(v_m==0)//dz^2
			d1=-3^0.5/3*(a*d+b*e-2*c*f)
			d2=3^0.5/6*(-a^2-b^2+2*c^2+d^2+e^2-2*f^2)
			d3=-3^0.5/3*(a*g+b*h-2*c*o)
			d4=-3^0.5/3*(d*g+e*h-2*f*o)
			d5=1/6*(a^2+b^2-2*c^2+d^2+e^2-2*(f^2+g^2+h^2-2*o^2))			
		elseif(v_m==-1)//dyz
			d1=c*e+b*f
			d2=b*c-e*f
			d3=c*h+b*o
			d4=f*h+e*o
			d5=-3^0.5/3*(b*c+e*f-2*h*o)

		elseif(v_m==-2)//dxy
			d1=b*d+a*c
			d2=a*b-d*e
			d3=b*g+a*h
			d4=e*g+d*h
			d5=-3^0.5/3*(a*b+d*e-2*g*h)
		endif
	///////////smple Rot. part at d-orb
	
	if(kind==1)
		realsp=abs(d1*PrepareRealSpH(v_l,-2,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+d2*PrepareRealSpH(v_l,2,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+d3*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+d4*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+d5*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					)
	elseif(kind==2)
		realsp=abs(d1*PrepareRealSpH(v_l,-2,kf,x,y,check_angle_or_k)+d2*PrepareRealSpH(v_l,2,kf,x,y,check_angle_or_k)\
					+d3*PrepareRealSpH(v_l,1,kf,x,y,check_angle_or_k)+d4*PrepareRealSpH(v_l,-1,kf,x,y,check_angle_or_k)\
					+d5*PrepareRealSpH(v_l,0,kf,x,y,check_angle_or_k)\
					)
	elseif(kind==3)
		realsp=abs(d1*PrepareRealSpH(v_l,-2,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+d2*PrepareRealSpH(v_l,2,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+d3*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+d4*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+d5*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					)
	
	else
		abort
	endif
	
	
/////f-orbital part//////////////
	elseif(v_l==3 && -abs(v_l)<=v_m && v_m<=v_l)
		variable/D f1, f2, f3, f4, f5, f6, f7
		if (v_m==3)//f_y(3x^2-y^2)
			f1 = 0.25*(3*b*a^2-6*d*e*a- b*(b^2+3*d^2-3*e^2) )//m=-3
			f2 = 0.25*(e^3+3*a^2*e-3*b^2*e-3*d^2*e+6*a*b*d)//m=3
			f3 = 6^0.5/2*(b*d*g+a*e*g+a*d*h-b*e*h)//m=2
			f4 = 6^0.5/4*(h*a^2+2*b*g*a-2*d*e*g-d^2*h +(e^2-b^2)*h )//m=-2
			f5 = -15^0.5/20*(3*b*a^2 +2*(d*e-4*g*h)*a -b*(b^2-d^2+e^2+4*g^2-4*h^2) )//m=-1
			f6 = 15^0.5/20*(-e*a^2-2*b*d*a-4*e*h^2 +e*(b^2-3*d^2+e^2+4*g^2) +8*d*g*h)//m=1
			f7 = -10^0.5/20*(2*h^3+3*a^2*h+3*d^2*h -3*(b^2+e^2+2*g^2)*h +6*a*b*g+6*d*e*g)//m=0
			
		elseif(v_m==2)//xyz
			f1 = -6^0.5/2*(-a*b*c+d*e*c+b*d*f+a*e*f)//m=-3
			f2 = 6^0.5/2*(b*c*d-e*f*d+a*c*e+a*b*f)//m=3
			f3 = c*e*g*+b*f*g+c*d*h+a*f*h+b*d*o+a*e*o//m=2
			f4 = b*c*g-e*f*g+a*c*h-d*f*h+a*b*o-d*e*o//m=-2
			f5 = -10^0.5/10*(c*d*e+b*d*f-4*c*g*h-4*b*g*o +a*(3*b*c+e*f-4*h*o) )//m=-1
			
			f6 = -10^0.5/10*(b*c*d+3*e*f*d-4*h*o*d+a*c*e +a*b*f-4*f*g*h-4*e*g*o)//m=1
			
			f7 = -15^0.5/5*(b*c*g+e*f*g-2*h*o*g+a*c*h+d*f*h+a*b*o+d*e*o)//m=0

		elseif(v_m==1)//f_yz^2
			f1 = -15^0.5/20*(b^3+a^2*b -(4*c^2+d^2+3*e^2-4*f^2)*b -2*a*d*e+8*c*e*f )//m=-3
			f2 = -15^0.5/20*(e*a^2+2*b*d*a+3*b^2*e-8*b*c*f -e*(4*c^2+d^2+e^2-4*f^2) )//m=3
			f3 = -10^0.5/10*(b*d*g+a*e*g+a*d*h+3*b*e*h-4*c*f*h-4*(c*e+b*f)*o)//m=2
			f4 = -10^0.5/20*(h*a^2+2*b*g*a-2*d*e*g+3*b^2*h-4*c^2*h-d^2*h-3*e^2*h+4*f^2*h-8*b*c*o+8*e*f*o)//m=-2
			f5 = 1/20*(3*b^3+3*a^2*b+(-12*c^2+d^2+3*e^2-4*f^2-4*g^2-12*h^2+16*o^2)*b+2*a*(d*e-4*g*h)-8*c*(e*f-4*h*o) )//m=-1
			f6 = 1/20*(3*e^3+a^2*e+3*b^2*e-4*c^2*e+3*d^2*e-12*f^2*e-4*g^2*e-12*h^2*e+16*o^2*e+2*a*b*d-8*b*c*f-8*d*g*h+32*f*h*o)//m=1
			
			f7 = 6^0.5/20*(-2*h^3+a^2*h+3*b^2*h-4*c^2*h+d^2*h+3*e^2*h-4*f^2*h-2*g^2*h+8*o^2*h+2*a*b*g+2*d*e*g-8*b*c*o-8*e*f*o)//m=0
		
		elseif(v_m==0)//f_z^3
			f1 = 10^0.5/20*(-3*c*a^2+6*d*f*a-3*b^2*c+6*b*e*f +c*(2*c^2+3*d^2+3*e^2-6*f^2) )//m=-3
			f2 = -10^0.5/20*(2*f^3+3*a^2*f+3*b^2*f-3*(2*c^2+d^2+e^2)*f+6*a*c*d+6*b*c*e )//m=3
			f3 = -15^0.5/5*(a*f*g+b*f*h+a*d*o+b*e*o+c*(d*g+e*h-2*f*o) )//m=2
			f4 = -15^0.5/10*( 2*(a*c*g+b*c*h-f*(d*g+e*h)) + (a^2+b^2-2*c^2-d^2-e^2+2*f^2)*o )//m=-2
			f5 = 6^0.5/20*(3*c*a^2+2*(d*f-4*g*o)*a +3*b^2*c+2*b*(e*f-4*h*o) +c*(-2*c^2+d^2+e^2-2*f^2-4*g^2-4*h^2+8*o^2) )//m=-1
			f6 = 6^0.5/20*(-2*f^3+a^2*f+b^2*f-2*c^2*f+3*d^2*f-4*g^2*f-4*h^2*f+8*o^2*f+2*a*c*d+2*b*c*e-8*d*g*o-8*e*h*o)//m=1
			f7 = 1/10*(4*o^3+3*(a^2+b^2-2*c^2+d^2+e^2-2*(f^2+g^2+h^2) )*o +6*(a*c*g+d*f*g+b*c*h+e*f*h) ) //m=0

		
		elseif(v_m==-1)//f_xz^2
			f1 = -15^0.5/20*(a^3 +(b^2-4*c^2-3*d^2-e^2+4*f^2)*a-2*b*d*e+8*c*d*f )//m=-3
			f2 = -15^0.5/20*(3*d*a^2+2*(b*e-4*c*f)*a-d*(-b^2+4*c^2+d^2+e^2-4*f^2) )//m=3
			f3 = -10^0.5/10*(3*a*d*g+b*e*g-4*c*f*g+b*d*h+a*e*h-4*(c*d+a*f)*o )//m=2
			f4 = -10^0.5/20*(3*g*a^2+2*b*h*a-8*c*o*a+b^2*g-4*c^2*g-3*d^2*g-e^2*g+4*f^2*g-2*d*e*h+8*d*f*o)//m=-2
			f5 = 1/20*(3*a^3 +(3*b^2-12*c^2+3*d^2+e^2-4*f^2-12*g^2-4*h^2+16*o^2)*a +2*b*d*e-8*c*d*f-8*b*g*h+32*c*g*o )//m=-1
			f6 = 1/20*(3*d^3+3*a^2*d+b^2*d-4*c^2*d+3*e^2*d-12*f^2*d-12*g^2*d-4*h^2*d+16*o^2*d+2*a*b*e-8*a*c*f-8*e*g*h+32*f*g*o)//m=1
			f7 = 6^0.5/20*(-2*g^3+3*a^2*g+b^2*g-4*c^2*g+3*d^2*g+e^2*g-4*f^2*g-2*h^2*g+8*o^2*g+2*a*b*h+2*d*e*h-8*a*c*o-8*d*f*o ) //m=0

		
		elseif(v_m==-2)//z(x^2-y^2)
			f1 = 6^0.5/4*(c*(a^2-b^2-d^2+e^2)-2*a*d*f+2*b*e*f )//m=-3
			f2 = 6^0.5/4*(f*a^2+2*c*d*a-2*b*c*e-b^2*f+(e^2-d^2)*f )//m=3
			f3 = c*d*g+a*f*g-c*e*h-b*f*h+a*d*o-b*e*o //m=2
			f4 = 0.5*(o*a^2+2*c*g*a-2*d*f*g-2*b*c*h+2*e*f*h-b^2*o-d^2*o+e^2*o)//m=-2
			f5 = -10^0.5/20*(3*c*a^2 +2*(d*f-4*g*o)*a -3*b^2*c +c*(d^2-e^2-4*g^2+4*h^2) +b*(8*h*o-2*e*f) )//m=-1
			f6 = -10^0.5/20*(f*a^2+2*c*d*a-4*f*g^2+4*f*h^2-2*b*c*e-b^2*f+3*d^2*f-3*e^2*f-8*d*g*o+8*e*h*o)//m=1
			f7 = 15^0.5/10*(o*a^2+2*c*g*a+2*d*f*g -2*(b*c+e*f)*h +d^2*o -(b^2+e^2+2*g^2-2*h^2)*o ) //m=0

		elseif(v_m==-3)//f_x(x^2-3y^2)ここ
			f1 = 0.25*(a^3 -3*(b^2+d^2-e^2)*a +6*b*d*e )//m=-3	
			f2 = 0.25*(3*d*a^2-6*b*e*a -d*(3*b^2+d^2-3*e^2) )//m=3
			f3 = -6^0.5/2*(-a*d*g+b*e*g+b*d*h+a*e*h) //m=2
			f4 = 6^0.5/4*( (a^2-b^2-d^2+e^2)*g -2*a*b*h+2*d*e*h)//m=-2
			f5 = -15^0.5/20*(a^3 +(-3*b^2+d^2-e^2-4*g^2+4*h^2)*a-2*b*d*e+8*b*g*h )//m=-1
			f6 = -15^0.5/20*(d^3+a^2*d-b^2*d-3*e^2*d-4*g^2*d+4*h^2*d-2*a*b*e+8*e*g*h)//m=1
			
			f7 = -10^0.5/20*(-2*g^3+3*a^2*g-3*b^2*g-3*e^2*g+6*h^2*g-6*a*b*h-6*d*e*h ) //m=0
		
		endif
		
		if(kind==1)
			realsp=abs(f1*PrepareRealSpH(v_l,-3,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+f2*PrepareRealSpH(v_l,3,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+f3*PrepareRealSpH(v_l,2,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+f4*PrepareRealSpH(v_l,-2,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+f5*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+f6*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+f7*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					)
		elseif(kind==2)
			realsp=abs(f1*PrepareRealSpH(v_l,-3,kf,x,y,check_angle_or_k)+f2*PrepareRealSpH(v_l,3,kf,x,y,check_angle_or_k)\
					+f3*PrepareRealSpH(v_l,2,kf,x,y,check_angle_or_k)+f4*PrepareRealSpH(v_l,-2,kf,x,y,check_angle_or_k)\
					+f5*PrepareRealSpH(v_l,-1,kf,x,y,check_angle_or_k)+f6*PrepareRealSpH(v_l,1,kf,x,y,check_angle_or_k)\
					+f7*PrepareRealSpH(v_l,0,kf,x,y,check_angle_or_k)\
					)
	
		elseif(kind==3)
			realsp=abs(f1*PrepareRealSpH(v_l,-3,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+f2*PrepareRealSpH(v_l,3,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+f3*PrepareRealSpH(v_l,2,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+f4*PrepareRealSpH(v_l,-2,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+f5*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)+f6*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					+f7*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,z,check_angle_or_k)\
					)
		else
			abort
		endif
		
		
	endif
	return realsp
end

////////////prepare/////////
function PrepareRealSpH(l,m,r,x,y,check_angle_or_k)
	variable l,m,r,x,y,check_angle_or_k
	variable Realspherical_harmonics
	 if(check_angle_or_k==0)//angle(DA30mode) if y=0 tiltmap mode
	 
	 	if(m<=-1)
	 		//Realspherical_harmonics= (-imag( (-1)^(-m) \
	 										//*sphericalHarmonics(l, abs(m), asin((sin(x*pi/180)^2+sin(y*pi/180)^2)^0.5), Polaris_nan2zero(sign(y*sin(((x*pi/180)^2+(y*pi/180)^2)^0.5))*atan(y/x) ) )\
	 										//-sphericalHarmonics(l, -abs(m), asin((sin(x*pi/180)^2+sin(y*pi/180)^2)^0.5), Polaris_nan2zero(sign(y*sin(((x*pi/180)^2+(y*pi/180)^2)^0.5))*atan(y/x) ) ) \
	 										//)/2^0.5)
	 			 		Realspherical_harmonics= (-imag( (-1)^(-m) \
	 										*sphericalHarmonics(l, abs(m), ((x*pi/180)^2+(y*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(y*pi/180))*acos(x/(x^2+y^2)^0.5*sin((x^2+y^2)^0.5*pi/180)/(x^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2+y^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2)^0.5)) )\
	 										-sphericalHarmonics(l, -abs(m),((x*pi/180)^2+(y*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(y*pi/180))*acos(x/(x^2+y^2)^0.5*sin((x^2+y^2)^0.5*pi/180)/(x^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2+y^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2)^0.5)) )  \
	 										)/(2^0.5))  
	 	elseif(m==0)
	 		//(1-x^-y^2)^0.5
	 		//realsp=abs(real(sphericalHarmonics(l, m, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))))*realsp
	 		Realspherical_harmonics=(real(sphericalHarmonics(l, m, ((x*pi/180)^2+(y*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(y*pi/180))*acos(x/(x^2+y^2)^0.5*sin((x^2+y^2)^0.5*pi/180)/(x^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2+y^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2)^0.5)) ) ))
	 	elseif(m>=1)
	 	
	 		Realspherical_harmonics=(real( (-1)^(m) \
	 										*sphericalHarmonics(l, abs(m), ((x*pi/180)^2+(y*pi/180)^2)^0.5,  Polaris_nan2zero(sign(sin(y*pi/180))*acos(x/(x^2+y^2)^0.5*sin((x^2+y^2)^0.5*pi/180)/(x^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2+y^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2)^0.5)) )\
	 										+sphericalHarmonics(l, -abs(m), ((x*pi/180)^2+(y*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(y*pi/180))*acos(x/(x^2+y^2)^0.5*sin((x^2+y^2)^0.5*pi/180)/(x^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2+y^2/(x^2+y^2)*sin((x^2+y^2)^0.5*pi/180)^2)^0.5)) ) \
	 										)/(2^0.5))
	 	endif
	 	
	 	
	 	
	 	
	 	
	 	
	 	
	 	
	 	
	 	
	 elseif(check_angle_or_k==1)//mometum(DA30 mode)
	 	if(m<=-1)
	 		Realspherical_harmonics= (-imag( (-1)^(-m) *sphericalHarmonics(l, abs(m), asin((x^2+y^2)^0.5/r), Polaris_nan2zero(sign(y)*acos(x/(x^2+y^2)^0.5)) )\
	 		-sphericalHarmonics(l, -abs(m), asin((x^2+y^2)^0.5/r), Polaris_nan2zero(sign(y)*acos(x/(x^2+y^2)^0.5))) )/(2^0.5)) 
	 	elseif(m==0)
	 		//(1-x^-y^2)^0.5
	 		//realsp=abs(real(sphericalHarmonics(l, m, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))))*realsp
	 		Realspherical_harmonics=(real(sphericalHarmonics(l, m, asin((x^2+y^2)^0.5/r), Polaris_nan2zero(sign(y)*acos(x/(x^2+y^2)^0.5)))))
	 	elseif(m>=1)
	 		Realspherical_harmonics=(real( (-1)^(m) *sphericalHarmonics(l, abs(m), asin((x^2+y^2)^0.5/r), Polaris_nan2zero(sign(y)*acos(x/(x^2+y^2)^0.5)))\
	 		+sphericalHarmonics(l, -abs(m), asin((x^2+y^2)^0.5/r), Polaris_nan2zero(sign(y)*acos(x/(x^2+y^2)^0.5))) )/(2^0.5))
		endif








	elseif(check_angle_or_k==2) //k tiltmapmode->破棄
		variable theta_x=asin(x/r)
		variable theta_y=0//atan(y/(r^2-x^2-y^2)^0.5)
	 	if(m<=-1)
	 		Realspherical_harmonics= (-imag( (-1)^(-m) *sphericalHarmonics(l, abs(m), ((theta_x)^2+(theta_y)^2)^0.5,0)\
	 		-sphericalHarmonics(l, -abs(m), ((theta_x)^2+(theta_y)^2)^0.5, 0) )/(2^0.5)) 
	 	elseif(m==0)
	 		//(1-x^-y^2)^0.5
	 		//realsp=abs(real(sphericalHarmonics(l, m, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))))*realsp
	 		Realspherical_harmonics=(real(sphericalHarmonics(l, m, ((theta_x)^2+(theta_y)^2)^0.5, 0)))
	 	elseif(m>=1)
	 		Realspherical_harmonics=(real( (-1)^(m) *sphericalHarmonics(l, abs(m), ((theta_x)^2+(theta_y)^2)^0.5, 0 )\
	 		+sphericalHarmonics(l, -abs(m), ((theta_x)^2+(theta_y)^2)^0.5, 0 ))/(2^0.5))
		endif
	 	
	 	
	 	
	 	
	 	
	 	
	 ///破棄
	elseif(check_angle_or_k==3) //k tiltmapmode(only tipe2)
		wave mw_y=root:pol_ARPES:misc:rot_matrix:rot_mat_y
		variable/D si_pol=mw_y[0][2]
		variable/D co_pol=mw_y[2][2]
	
		wave mw_z=root:pol_ARPES:misc:rot_matrix:rot_mat_z
		variable/D si_azi=mw_z[1][0]
		variable/D co_azi=mw_z[1][1]
		
		nvar v_tilt=root:pol_ARPES:misc:v_tilt
		
		//variable alpha=180/pi*asin( (si_pol*sqrt(r^2-(si_azi*x-co_azi*y)^2) +co_pol*(si_azi*x-co_azi*y))/r )
		//variable bet=ti_ofs+180/pi*atan((co_azi*x+si_azi*y)/(r^2-x^2-y^2)^0.5 )
		variable/D alpha=-180/pi*asin( (co_azi*co_pol-si_azi*si_pol*(co_azi*y+si_azi*x)/(r^2-x^2-y^2+(co_azi*y+si_azi*x)^2)^0.5)*x/r\
									+(si_azi*co_pol+co_azi*si_pol*(co_azi*y+si_azi*x)/(r^2-x^2-y^2+(co_azi*y+si_azi*x)^2)^0.5)*y/r\
									-si_pol*(r^2-x^2-y^2)/(r^2-x^2-y^2+(co_azi*y+si_azi*x)^2)^0.5/r)
		variable/D bet=v_tilt+180/pi*atan((co_azi*y+si_azi*x)/(r^2-x^2-y^2)^0.5 )
		if(m<0)
	 	 	Realspherical_harmonics=(-imag( (-1)^(-m) *sphericalHarmonics(l, abs(m),((alpha*pi/180)^2+(bet*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(bet*pi/180))*acos(sin(alpha*pi/180)/(sin(alpha*pi/180)^2+sin(bet*pi/180)^2)^0.5)) )\
	 						-sphericalHarmonics(l, -abs(m), ((alpha*pi/180)^2+(bet*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(bet*pi/180))*acos(sin(alpha*pi/180)/(sin(alpha*pi/180)^2+sin(bet*pi/180)^2)^0.5)) ))/2^0.5)
	 	//Realspherical_harmonics= (-imag( (-1)^(-m) *sphericalHarmonics(l, abs(m), asin((y^2+0^2)^0.5/x), sign(r)*acos(y/(y^2+r^2)^0.5))-sphericalHarmonics(l, -abs(m), asin((y^2+r^2)^0.5/x), sign(r)*acos(y/(y^2+r^2)^0.5)) )/2^0.5) 
	 	elseif(m==0)
	 	//(1-x^-y^2)^0.5
	 	//realsp=abs(real(sphericalHarmonics(l, m, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))))*realsp//sign(cos(x*pi/180)*sin(y*pi/180))*
	 	//Realspherical_harmonics_test=(real(sphericalHarmonics(l, m, acos(cos(x*pi/180)*cos(y*pi/180)), asin((cos(x*pi/180)*sin(y*pi/180))/((cos(x*pi/180))^2*(sin(y*pi/180))^2+(sin(x*pi/180))^2)^0.5  )) ))
	 		Realspherical_harmonics=(real(sphericalHarmonics(l, m,((alpha*pi/180)^2+(bet*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(bet*pi/180))*acos(sin(alpha*pi/180)/(sin(alpha*pi/180)^2+sin(bet*pi/180)^2)^0.5)) ) ))

	 		
	 	elseif(m>0)
	 	
	 		Realspherical_harmonics=(real( (-1)^(-m) *sphericalHarmonics(l, abs(m),((alpha*pi/180)^2+(bet*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(bet*pi/180))*acos(sin(alpha*pi/180)/(sin(alpha*pi/180)^2+sin(bet*pi/180)^2)^0.5)) )\
	 						+sphericalHarmonics(l, -abs(m), ((alpha*pi/180)^2+(bet*pi/180)^2)^0.5, Polaris_nan2zero(sign(sin(bet*pi/180))*acos(sin(alpha*pi/180)/(sin(alpha*pi/180)^2+sin(bet*pi/180)^2)^0.5)) ))/2^0.5)
	 	endif
	 
	endif	
	return Realspherical_harmonics
end

//////////////
//function PrepareRealSpH_tilt_ang(l,m,x,y,r)
//	variable l,m,x,y,r

//	variable Realspherical_harmonics
	 
//angle tiltmapmode(only tipe2)
//		 if(m<0)
//	 	 	Realspherical_harmonics=(-imag( (-1)^(-m) *\
//	 	 									sphericalHarmonics(l, abs(m), acos(cos(x*pi/180)), 0)\
//	 										-sphericalHarmonics(l, -abs(m), acos(cos(x*pi/180)*cos(y*pi/180)), 0) )/2^0.5)
	 	//Realspherical_harmonics= (-imag( (-1)^(-m) *sphericalHarmonics(l, abs(m), asin((y^2+0^2)^0.5/x), sign(r)*acos(y/(y^2+r^2)^0.5))-sphericalHarmonics(l, -abs(m), asin((y^2+r^2)^0.5/x), sign(r)*acos(y/(y^2+r^2)^0.5)) )/2^0.5) 
//	 	elseif(m==0)
	 	//(1-x^-y^2)^0.5
	 	//realsp=abs(real(sphericalHarmonics(l, m, acos(z/(x^2+y^2+z^2)^0.5), sign(y)*acos(x/(x^2+y^2)^0.5))))*realsp//sign(cos(x*pi/180)*sin(y*pi/180))*
	 	//Realspherical_harmonics_test=(real(sphericalHarmonics(l, m, acos(cos(x*pi/180)*cos(y*pi/180)), asin((cos(x*pi/180)*sin(y*pi/180))/((cos(x*pi/180))^2*(sin(y*pi/180))^2+(sin(x*pi/180))^2)^0.5  )) ))
//	 		Realspherical_harmonics=(real(sphericalHarmonics(l, m,acos(cos(x*pi/180)*cos(y*pi/180)) , 0) ))
//	 	elseif(m>0)
	 	
//	 		Realspherical_harmonics=(real( (-1)^(-m) *sphericalHarmonics(l, abs(m), acos(cos(x*pi/180)*cos(y*pi/180)), 0)\
//	 						+sphericalHarmonics(l, -abs(m), acos(cos(x*pi/180)*cos(y*pi/180)), 0))/2^0.5)
//	 	endif
//	 	return Realspherical_harmonics
//end


///wave作成用関数
function PrepareWave(kind, pointE, Est, Edel, pointkx, kx_st, kx_del,pointky,ky_st, ky_del)
	 variable kind, pointE, Est, Edel, pointkx,  kx_st, kx_del,pointky,ky_st, ky_del
	 
	
	if(kind==1)///E-kの時(検出面)
		Make/D/N=(pointE,pointkx)/D/O pre_w ///x=E,y=kx (inomacro仕様)
		wave pre_w
		setscale/p x, Est, Edel,pre_w ////E
		setscale/p y, kx_st, kx_del,pre_w ///kx
	elseif(kind==2)//k-k(フェルミ面のみ)
		Make/D/N=(pointkx,pointky)/D/O pre_w
		wave pre_w
		setscale/p x, kx_st, kx_del,pre_w ////kx
		setscale/p y, ky_st, ky_del,pre_w ///ky
	elseif(kind==3)//E(kx,ky)volデータ//DA30 mode
		Make/D/N=(pointE,pointkx,pointky)/D/O pre_w
		wave pre_w
		setscale/p x, Est, Edel,pre_w ////E
		setscale/p y, kx_st, kx_del,pre_w ///kx
		setscale/p z, ky_st, ky_del,pre_w ///ky
	else
	abort "please select which a kind ARPES sim. deta"
	endif
end  


Function Polaris_Wigner_Rot_Matrix_tiltmap(v_l,v_m,v_hn,v_W, v_pol,v_tilt,v_azi,check_angle_or_k, kind, x,y,z)///実装用
	variable/D v_l,v_m,v_hn,v_W, v_pol,v_azi,kind,check_angle_or_k,v_tilt, x,y,z //deg
	variable realsp
	v_tilt=v_tilt*(-1)
	
	//variable check_angle_or_k=1
	
	variable E_b
	variable/D kf=0.512*(v_hn-v_W+E)^0.5//(エネルギー方向の変数を入れる)
	variable/D a, b, c, d, e, f, g, h, o //回転行列の要素
	string l_s, m_s, sph
	
	if(kind==1 && check_angle_or_k==0)
 		a= Cos(v_azi*pi/180)*Cos(v_pol*pi/180) - Sin(v_azi*pi/180)*Sin(v_pol*pi/180)*Sin((v_tilt)*pi/180)
		b= -Cos((v_tilt)*pi/180)* Sin(v_azi*pi/180)
		c= Cos(v_azi*pi/180)* Sin(v_pol*pi/180) + Cos(v_pol*pi/180)* Sin(v_azi*pi/180)* Sin((v_tilt)*pi/180)
		d= Cos(v_pol*pi/180)* Sin(v_azi*pi/180) + Cos(v_azi*pi/180)* Sin(v_pol*pi/180)* Sin((v_tilt)*pi/180)
	
		e= Cos(v_azi*pi/180)* Cos((v_tilt)*pi/180) 
		f= Sin(v_azi*pi/180)* Sin(v_pol*pi/180) - Cos(v_azi*pi/180)* Cos(v_pol*pi/180)* Sin((v_tilt)*pi/180)
		g= -Cos((v_tilt)*pi/180)* Sin(v_pol*pi/180)
		h= Sin((v_tilt)*pi/180)
		o= Cos(v_pol*pi/180)* Cos((v_tilt)*pi/180)

	elseif(kind==2 && check_angle_or_k==0)
 		a= Cos(v_azi*pi/180)*Cos(v_pol*pi/180) - Sin(v_azi*pi/180)*Sin(v_pol*pi/180)*Sin(-y*pi/180)
		b= -Cos(-y*pi/180)* Sin(v_azi*pi/180)
		c= Cos(v_azi*pi/180)* Sin(v_pol*pi/180) + Cos(v_pol*pi/180)* Sin(v_azi*pi/180)* Sin(-y*pi/180)
		d= Cos(v_pol*pi/180)* Sin(v_azi*pi/180) + Cos(v_azi*pi/180)* Sin(v_pol*pi/180)* Sin(-y*pi/180)
	
		e= Cos(v_azi*pi/180)* Cos(-y*pi/180) 
		f= Sin(v_azi*pi/180)* Sin(v_pol*pi/180) - Cos(v_azi*pi/180)* Cos(v_pol*pi/180)* Sin(-y*pi/180)
		g= -Cos(-y*pi/180)* Sin(v_pol*pi/180)
		h= Sin(-y*pi/180)
		o= Cos(v_pol*pi/180)* Cos(-y*pi/180)
		
	elseif(kind==3 && check_angle_or_k==0)
	   a= Cos(v_azi*pi/180)*Cos(v_pol*pi/180) - Sin(v_azi*pi/180)*Sin(v_pol*pi/180)*Sin(-z*pi/180)
		b= -Cos(y*pi/180)* Sin(v_azi*pi/180)
		c= Cos(v_azi*pi/180)* Sin(v_pol*pi/180) + Cos(v_pol*pi/180)* Sin(v_azi*pi/180)* Sin(-z*pi/180)
		d= Cos(v_pol*pi/180)* Sin(v_azi*pi/180) + Cos(v_azi*pi/180)* Sin(v_pol*pi/180)* Sin(-z*pi/180)
	
		e= Cos(v_azi*pi/180)* Cos(y*pi/180) 
		f= Sin(v_azi*pi/180)* Sin(v_pol*pi/180) - Cos(v_azi*pi/180)* Cos(v_pol*pi/180)* Sin(-z*pi/180)
		g= -Cos(y*pi/180)* Sin(v_pol*pi/180)
		h= Sin(y*pi/180)
		o= Cos(v_pol*pi/180)* Cos(-z*pi/180)
		
	elseif(kind==1 && check_angle_or_k==1)
	
	//E-kx
	E_b=x
		a= Cos(v_azi*pi/180)*Cos(v_pol*pi/180) - Sin(v_azi*pi/180)*Sin(v_pol*pi/180)*v_tilt/(kf^2-y^2)^0.5
		b= -(kf^2-y^2-v_tilt^2)^0.5/(kf^2-y^2)^0.5* Sin(v_azi*pi/180)
		c= Cos(v_azi*pi/180)* Sin(v_pol*pi/180) + Cos(v_pol*pi/180)* Sin(v_azi*pi/180)* v_tilt/(kf^2-y^2)^0.5
		d= Cos(v_pol*pi/180)* Sin(v_azi*pi/180) + Cos(v_azi*pi/180)* Sin(v_pol*pi/180)* v_tilt/(kf^2-y^2)^0.5
	
		e= Cos(v_azi*pi/180)* (kf^2-y^2-v_tilt^2)^0.5/(kf^2-y^2)^0.5 
		f= Sin(v_azi*pi/180)* Sin(v_pol*pi/180) - Cos(v_azi*pi/180)* Cos(v_pol*pi/180)*v_tilt/(kf^2-y^2)^0.5
		g= -(kf^2-y^2-v_tilt^2)^0.5/(kf^2-y^2)^0.5* Sin(v_pol*pi/180)
		h= v_tilt/(kf^2-y^2)^0.5
		o= Cos(v_pol*pi/180)* (kf^2-y^2-v_tilt^2)^0.5/(kf^2-y^2)^0.5
		
		
		
	elseif(kind==2 && check_angle_or_k==1)
	
	//constmapのみ
		a= Cos(v_azi*pi/180)*Cos(v_pol*pi/180) - Sin(v_azi*pi/180)*Sin(v_pol*pi/180)*(-y)/(kf^2-x^2)^0.5//*sign(sin(y*pi/180))
		b= -(kf^2-x^2-y^2)^0.5/(kf^2-x^2)^0.5* Sin(v_azi*pi/180)//*sign(x)
		c= Cos(v_azi*pi/180)* Sin(v_pol*pi/180) + Cos(v_pol*pi/180)* Sin(v_azi*pi/180)* (-y)/(kf^2-x^2)^0.5//*sign(sin(y*pi/180))
		d= Cos(v_pol*pi/180)* Sin(v_azi*pi/180) + Cos(v_azi*pi/180)* Sin(v_pol*pi/180)* (-y)/(kf^2-x^2)^0.5//*sign(sin(y*pi/180))
	
		e= Cos(v_azi*pi/180)* (kf^2-x^2-y^2)^0.5/(kf^2-x^2)^0.5 //*sign(x)
		f= Sin(v_azi*pi/180)* Sin(v_pol*pi/180) - Cos(v_azi*pi/180)* Cos(v_pol*pi/180)//*y/(kf^2-x^2)^0.5*sign(sin(y*pi/180))
		g= -(kf^2-x^2-y^2)^0.5/(kf^2-x^2)^0.5* Sin(v_pol*pi/180)//*sign(x)
		h= (-y)/(kf^2-x^2)^0.5//*sign(sin(y*pi/180))
		o= Cos(v_pol*pi/180)* (kf^2-x^2-y^2)^0.5/(kf^2-x^2)^0.5//*sign(x)
		
	elseif(kind==3 && check_angle_or_k==1)//vol
		E_b=x
		a= Cos(v_azi*pi/180)*Cos(v_pol*pi/180) - Sin(v_azi*pi/180)*Sin(v_pol*pi/180)*(-z)/(kf^2-y^2)^0.5
		b= -(kf^2-y^2-z^2)^0.5/(kf^2-y^2)^0.5* Sin(v_azi*pi/180)
		c= Cos(v_azi*pi/180)* Sin(v_pol*pi/180) + Cos(v_pol*pi/180)* Sin(v_azi*pi/180)* (-z)/(kf^2-y^2)^0.5
		d= Cos(v_pol*pi/180)* Sin(v_azi*pi/180) + Cos(v_azi*pi/180)* Sin(v_pol*pi/180)* (-z)/(kf^2-y^2)^0.5
	
		e= Cos(v_azi*pi/180)* (kf^2-y^2-z^2)^0.5/(kf^2-y^2)^0.5 
		f= Sin(v_azi*pi/180)* Sin(v_pol*pi/180) - Cos(v_azi*pi/180)* Cos(v_pol*pi/180)*(-z)/(kf^2-y^2)^0.5
		g= -(kf^2-y^2-z^2)^0.5/(kf^2-y^2)^0.5* Sin(v_pol*pi/180)
		h= (-z)/(kf^2-y^2)^0.5
		o= Cos(v_pol*pi/180)* (kf^2-y^2-z^2)^0.5/(kf^2-y^2)^0.5
	endif
	/////////////関数代入パート///////////////////
	if(v_l==0 && v_m==0) //s-orb nochange
		if(kind==1)
			realsp=abs(PrepareRealSpH(0,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k) )
		elseif(kind==2 && check_angle_or_k==0)
			realsp=abs(PrepareRealSpH(0,0,kf,x,0,check_angle_or_k) )
		elseif(kind==3)
			realsp=abs(PrepareRealSpH(0,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k) )
		else
			Abort
		endif
	elseif(v_l==1 && -abs(v_l)<=v_m && v_m<=v_l) ///p-orb part
     //////px,py, pz////選択パート
		variable/D t, u, v
		if(v_m==1)//px
			t=a; u=d; v=g
		elseif(v_m==-1)//py
			t=b; u=e; v=h
		else //pz
			t=c; u=f; v=o
		endif
	///////////////////sample Rot. part////////////
	///tilymap時サンプルに対する検出面が変わってしまう問題をどうするか? =>for分で回す(脳筋)
	//		解決方法1. tiltabgele をfor文を回して連続変化させ、その都度ky(y)=kf*cos(pi/180*polar)*sin(pi/180*tilt)を代入:polar=constとする
	//ただし、データ点の間隔が一定で無くなる。=>kyを連続変化させれば良い設定すればよくない？
	
	//解決方法2. 計算して解析的に求める///むずい
		//realsp=abs(1/3*(t*(real( (-1)^(-1) *sphericalHarmonics(l, abs(1), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5))+sphericalHarmonics(l, -abs(1), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5)) )/2^0.5) \
				//+v*(real(sphericalHarmonics(l, 0, asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5)))) \
				//+u*(-imag( (-1)^(-(-1)) *sphericalHarmonics(l, abs(-1), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5))-sphericalHarmonics(l, -abs(-1), asin((x^2+y^2)^0.5/kf), sign(y)*acos(x/(x^2+y^2)^0.5)) )/2^0.5) )\
				//)
				
	////////////if文にてkx-kyかE-kかvolかtiltmapかを判断			
			
	switch(kind)
		case 1:	//E-kx
			realsp=abs((t*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+v*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+u*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)))
		break
			
		case 2://only FS
			realsp=abs((t*PrepareRealSpH(v_l,1,kf,x,0,check_angle_or_k)+v*PrepareRealSpH(v_l,0,kf,x,0,check_angle_or_k)+u*PrepareRealSpH(v_l,-1,kf,x,0,check_angle_or_k)))
			//realsp=abs(PrepareRealSpH(v_l,0,kf,x,y,check_angle_or_k))
			
			//realsp=abs(1/3*(t*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+v*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+u*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)))
		break
		case 3://3Dvol
			realsp=abs((t*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W-x)^0.5,y,0,check_angle_or_k)+v*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W-x)^0.5,y,0,check_angle_or_k)+u*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W-0)^0.5,y,0,check_angle_or_k)))
		break
	endswitch
		
		
////////d-orbパート//////////////////////////	
	elseif(v_l==2 && -abs(v_l)<=v_m && v_m<=v_l) ///d-orb part
		variable/D d1, d2, d3, d4, d5
		if(v_m==2)//dx^2-y^2 rot
			d1=a*d-b*e //dxy'
			d2=0.5*(a^2-b^2-d^2+e^2) //dx^2-y^2'
			d3=a*g-b*h //dzx'
			d4=d*g-e*h //dyz'
			d5=-3^0.5/6*(a^2-b^2+d^2-e^2-2*g^2+2*h^2) //dz^2'
			
		elseif(v_m==1)//dzx
			d1=c*d+a*f //dxy'
			d2=a*c-d*f
			d3=c*g+a*o
			d4=f*g+d*o
			d5=-3^0.5/3*(a*c+d*f-2*g*o)
		elseif(v_m==0)//dz^2
			d1=-3^0.5/3*(a*d+b*e-2*c*f)
			d2=3^0.5/6*(-a^2-b^2+2*c^2+d^2+e^2-2*f^2)
			d3=-3^0.5/3*(a*g+b*h-2*c*o)
			d4=-3^0.5/3*(d*g+e*h-2*f*o)
			d5=1/6*(a^2+b^2-2*c^2+d^2+e^2-2*(f^2+g^2+h^2-2*o^2))			
		elseif(v_m==-1)//dyz
			d1=c*e+b*f
			d2=b*c-e*f
			d3=c*h+b*o
			d4=f*h+e*o
			d5=-3^0.5/3*(b*c+e*f-2*h*o)

		elseif(v_m==-2)//dxy
			d1=b*d+a*c
			d2=a*b-d*e
			d3=b*g+a*h
			d4=e*g+d*h
			d5=-3^0.5/3*(a*b+d*e-2*g*h)
		endif
	///////////smple Rot. part at d-orb
	if(kind==1)
		realsp=abs(d1*PrepareRealSpH(v_l,-2,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+d2*PrepareRealSpH(v_l,2,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+d3*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+d4*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+d5*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					)
	elseif(kind==2)
		realsp=abs(d1*PrepareRealSpH(v_l,-2,kf,x,0,check_angle_or_k)+d2*PrepareRealSpH(v_l,2,kf,x,0,check_angle_or_k)\
					+d3*PrepareRealSpH(v_l,1,kf,x,0,check_angle_or_k)+d4*PrepareRealSpH(v_l,-1,kf,x,0,check_angle_or_k)\
					+d5*PrepareRealSpH(v_l,0,kf,x,0,check_angle_or_k)\
					)
	elseif(kind==3)
		realsp=abs(d1*PrepareRealSpH(v_l,-2,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+d2*PrepareRealSpH(v_l,2,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+d3*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+d4*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+d5*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					)
	
	else
		abort
	endif
/////f-orbital part//////////////
	elseif(v_l==3 && -abs(v_l)<=v_m && v_m<=v_l)
		variable/D f1, f2, f3, f4, f5, f6, f7
		if (v_m==3)//f_y(3x^2-y^2)
			f1 = 0.25*(3*b*a^2-6*d*e*a- b*(b^2+3*d^2-3*e^2) )//m=-3
			f2 = 0.25*(e^3+3*a^2*e-3*b^2*e-3*d^2*e+6*a*b*d)//m=3
			f3 = 6^0.5/2*(b*d*g+a*e*g+a*d*h-b*e*h)//m=2
			f4 = 6^0.5/4*(h*a^2+2*b*g*a-2*d*e*g-d^2*h +(e^2-b^2)*h )//m=-2
			f5 = -15^0.5/20*(3*b*a^2 +2*(d*e-4*g*h)*a -b*(b^2-d^2+e^2+4*g^2-4*h^2) )//m=-1
			f6 = 15^0.5/20*(-e*a^2-2*b*d*a-4*e*h^2 +e*(b^2-3*d^2+e^2+4*g^2) +8*d*g*h)//m=1
			f7 = -10^0.5/20*(2*h^3+3*a^2*h+3*d^2*h -3*(b^2+e^2+2*g^2)*h +6*a*b*g+6*d*e*g)//m=0
			
		elseif(v_m==2)//xyz
			f1 = -6^0.5/2*(-a*b*c+d*e*c+b*d*f+a*e*f)//m=-3
			f2 = 6^0.5/2*(b*c*d-e*f*d+a*c*e+a*b*f)//m=3
			f3 = c*e*g*+b*f*g+c*d*h+a*f*h+b*d*o+a*e*o//m=2
			f4 = b*c*g-e*f*g+a*c*h-d*f*h+a*b*o-d*e*o//m=-2
			f5 = -10^0.5/10*(c*d*e+b*d*f-4*c*g*h-4*b*g*o +a*(3*b*c+e*f-4*h*o) )//m=-1
			
			f6 = -10^0.5/10*(b*c*d+3*e*f*d-4*h*o*d+a*c*e +a*b*f-4*f*g*h-4*e*g*o)//m=1
			
			f7 = -15^0.5/5*(b*c*g+e*f*g-2*h*o*g+a*c*h+d*f*h+a*b*o+d*e*o)//m=0

		elseif(v_m==1)//f_yz^2
			f1 = -15^0.5/20*(b^3+a^2*b -(4*c^2+d^2+3*e^2-4*f^2)*b -2*a*d*e+8*c*e*f )//m=-3
			f2 = -15^0.5/20*(e*a^2+2*b*d*a+3*b^2*e-8*b*c*f -e*(4*c^2+d^2+e^2-4*f^2) )//m=3
			f3 = -10^0.5/10*(b*d*g+a*e*g+a*d*h+3*b*e*h-4*c*f*h-4*(c*e+b*f)*o)//m=2
			f4 = -10^0.5/20*(h*a^2+2*b*g*a-2*d*e*g+3*b^2*h-4*c^2*h-d^2*h-3*e^2*h+4*f^2*h-8*b*c*o+8*e*f*o)//m=-2
			f5 = 1/20*(3*b^3+3*a^2*b+(-12*c^2+d^2+3*e^2-4*f^2-4*g^2-12*h^2+16*o^2)*b+2*a*(d*e-4*g*h)-8*c*(e*f-4*h*o) )//m=-1
			f6 = 1/20*(3*e^3+a^2*e+3*b^2*e-4*c^2*e+3*d^2*e-12*f^2*e-4*g^2*e-12*h^2*e+16*o^2*e+2*a*b*d-8*b*c*f-8*d*g*h+32*f*h*o)//m=1
			
			f7 = 6^0.5/20*(-2*h^3+a^2*h+3*b^2*h-4*c^2*h+d^2*h+3*e^2*h-4*f^2*h-2*g^2*h+8*o^2*h+2*a*b*g+2*d*e*g-8*b*c*o-8*e*f*o)//m=0
		
		elseif(v_m==0)//f_z^3
			f1 = 10^0.5/20*(-3*c*a^2+6*d*f*a-3*b^2*c+6*b*e*f +c*(2*c^2+3*d^2+3*e^2-6*f^2) )//m=-3
			f2 = -10^0.5/20*(2*f^3+3*a^2*f+3*b^2*f-3*(2*c^2+d^2+e^2)*f+6*a*c*d+6*b*c*e )//m=3
			f3 = -15^0.5/5*(a*f*g+b*f*h+a*d*o+b*e*o+c*(d*g+e*h-2*f*o) )//m=2
			f4 = -15^0.5/10*( 2*(a*c*g+b*c*h-f*(d*g+e*h)) + (a^2+b^2-2*c^2-d^2-e^2+2*f^2)*o )//m=-2
			f5 = 6^0.5/20*(3*c*a^2+2*(d*f-4*g*o)*a +3*b^2*c+2*b*(e*f-4*h*o) +c*(-2*c^2+d^2+e^2-2*f^2-4*g^2-4*h^2+8*o^2) )//m=-1
			f6 = 6^0.5/20*(-2*f^3+a^2*f+b^2*f-2*c^2*f+3*d^2*f-4*g^2*f-4*h^2*f+8*o^2*f+2*a*c*d+2*b*c*e-8*d*g*o-8*e*h*o)//m=1
			f7 = 1/10*(4*o^3+3*(a^2+b^2-2*c^2+d^2+e^2-2*(f^2+g^2+h^2) )*o +6*(a*c*g+d*f*g+b*c*h+e*f*h) ) //m=0

		
		elseif(v_m==-1)//f_xz^2
			f1 = -15^0.5/20*(a^3 +(b^2-4*c^2-3*d^2-e^2+4*f^2)*a-2*b*d*e+8*c*d*f )//m=-3
			f2 = -15^0.5/20*(3*d*a^2+2*(b*e-4*c*f)*a-d*(-b^2+4*c^2+d^2+e^2-4*f^2) )//m=3
			f3 = -10^0.5/10*(3*a*d*g+b*e*g-4*c*f*g+b*d*h+a*e*h-4*(c*d+a*f)*o )//m=2
			f4 = -10^0.5/20*(3*g*a^2+2*b*h*a-8*c*o*a+b^2*g-4*c^2*g-3*d^2*g-e^2*g+4*f^2*g-2*d*e*h+8*d*f*o)//m=-2
			f5 = 1/20*(3*a^3 +(3*b^2-12*c^2+3*d^2+e^2-4*f^2-12*g^2-4*h^2+16*o^2)*a +2*b*d*e-8*c*d*f-8*b*g*h+32*c*g*o )//m=-1
			f6 = 1/20*(3*d^3+3*a^2*d+b^2*d-4*c^2*d+3*e^2*d-12*f^2*d-12*g^2*d-4*h^2*d+16*o^2*d+2*a*b*e-8*a*c*f-8*e*g*h+32*f*g*o)//m=1
			f7 = 6^0.5/20*(-2*g^3+3*a^2*g+b^2*g-4*c^2*g+3*d^2*g+e^2*g-4*f^2*g-2*h^2*g+8*o^2*g+2*a*b*h+2*d*e*h-8*a*c*o-8*d*f*o ) //m=0

		
		elseif(v_m==-2)//z(x^2-y^2)
			f1 = 6^0.5/4*(c*(a^2-b^2-d^2+e^2)-2*a*d*f+2*b*e*f )//m=-3
			f2 = 6^0.5/4*(f*a^2+2*c*d*a-2*b*c*e-b^2*f+(e^2-d^2)*f )//m=3
			f3 = c*d*g+a*f*g-c*e*h-b*f*h+a*d*o-b*e*o //m=2
			f4 = 0.5*(o*a^2+2*c*g*a-2*d*f*g-2*b*c*h+2*e*f*h-b^2*o-d^2*o+e^2*o)//m=-2
			f5 = -10^0.5/20*(3*c*a^2 +2*(d*f-4*g*o)*a -3*b^2*c +c*(d^2-e^2-4*g^2+4*h^2) +b*(8*h*o-2*e*f) )//m=-1
			f6 = -10^0.5/20*(f*a^2+2*c*d*a-4*f*g^2+4*f*h^2-2*b*c*e-b^2*f+3*d^2*f-3*e^2*f-8*d*g*o+8*e*h*o)//m=1
			f7 = 15^0.5/10*(o*a^2+2*c*g*a+2*d*f*g -2*(b*c+e*f)*h +d^2*o -(b^2+e^2+2*g^2-2*h^2)*o ) //m=0

		elseif(v_m==-3)//f_x(x^2-3y^2)ここ
			f1 = 0.25*(a^3 -3*(b^2+d^2-e^2)*a +6*b*d*e )//m=-3	
			f2 = 0.25*(3*d*a^2-6*b*e*a -d*(3*b^2+d^2-3*e^2) )//m=3
			f3 = -6^0.5/2*(-a*d*g+b*e*g+b*d*h+a*e*h) //m=2
			f4 = 6^0.5/4*( (a^2-b^2-d^2+e^2)*g -2*a*b*h+2*d*e*h)//m=-2
			f5 = -15^0.5/20*(a^3 +(-3*b^2+d^2-e^2-4*g^2+4*h^2)*a-2*b*d*e+8*b*g*h )//m=-1
			f6 = -15^0.5/20*(d^3+a^2*d-b^2*d-3*e^2*d-4*g^2*d+4*h^2*d-2*a*b*e+8*e*g*h)//m=1
			
			f7 = -10^0.5/20*(-2*g^3+3*a^2*g-3*b^2*g-3*e^2*g+6*h^2*g-6*a*b*h-6*d*e*h ) //m=0
		
		endif
		
		if(kind==1)
			realsp=abs(f1*PrepareRealSpH(v_l,-3,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+f2*PrepareRealSpH(v_l,3,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+f3*PrepareRealSpH(v_l,2,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+f4*PrepareRealSpH(v_l,-2,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+f5*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+f6*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+f7*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					)
		elseif(kind==2)
			realsp=abs(f1*PrepareRealSpH(v_l,-3,kf,x,y,check_angle_or_k)+f2*PrepareRealSpH(v_l,3,kf,x,0,check_angle_or_k)\
					+f3*PrepareRealSpH(v_l,2,kf,x,0,check_angle_or_k)+f4*PrepareRealSpH(v_l,-2,kf,x,0,check_angle_or_k)\
					+f5*PrepareRealSpH(v_l,-1,kf,x,0,check_angle_or_k)+f6*PrepareRealSpH(v_l,1,kf,x,0,check_angle_or_k)\
					+f7*PrepareRealSpH(v_l,0,kf,x,0,check_angle_or_k)\
					)
	
		elseif(kind==3)
			realsp=abs(f1*PrepareRealSpH(v_l,-3,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+f2*PrepareRealSpH(v_l,3,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+f3*PrepareRealSpH(v_l,2,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+f4*PrepareRealSpH(v_l,-2,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+f5*PrepareRealSpH(v_l,-1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)+f6*PrepareRealSpH(v_l,1,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					+f7*PrepareRealSpH(v_l,0,0.512*(v_hn-v_W+x)^0.5,y,0,check_angle_or_k)\
					)
		else
			abort
		endif
		
		
	endif
	return realsp
end







Function ButtonProc_makewave(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	SetDataFolder root:Polaris:output
	
	NVar v_pol_st=root:Polaris:misc:v_theta_s
	NVAR v_pol_del=root:Polaris:misc:v_polar_del
	NVAR v_pol_point=root:Polaris:misc:v_polar_point
	
	//NVAR v_pol_en=root:Polaris:misc:v_theta_e
	
	NVAR v_tilt_st=root:Polaris:misc:v_tilt_s
	NVAR v_tilt_del=root:Polaris:misc:v_tilt_del
	NVAR v_tilt_point=root:Polaris:misc:v_tilt_point
	
	//NVAR v_tilt_en=root:Polaris:misc:v_tilt_e
	
	NVAR v_Est=root:Polaris:misc:v_Est
	NVAR v_Edel=root:Polaris:misc:v_Edel
	NVAR v_Epoint=root:Polaris:misc:v_Epoint
	
	NVAR v_kind=root:Polaris:misc:v_kind
	
	//NVAR v_Een=root:Polaris:misc:v_Een


	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			break
		case -1: // control being killed
			break
	endswitch

	return PrepareWave(v_kind, v_Epoint,v_Est,v_Edel, v_pol_point, v_pol_st, v_pol_del,v_tilt_point, v_tilt_st, v_tilt_del)
End



Function Polaris_PopMenu_orbital_selection_Proc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	NVAR v_l=root:Polaris:misc:v_l
	NVAR v_s=root:Polaris:misc:v_s
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			if(popNum==1)
				v_l=0; v_s=0
			elseif(popNum==2)
				v_l=1;v_s=1
			elseif(popNum==3)
				v_l=1;v_s=-1
			elseif(popNum==4)
				v_l=1;v_s=0
			elseif(popNum==5)
				v_l=2;v_s=-2
			elseif(popNum==6)
				v_l=2;v_s=2
			elseif(popNum==7)
				v_l=2;v_s=1
			elseif(popNum==8)
				v_l=2;v_s=-1
			elseif(popNum==9)
				v_l=2;v_s=0
			elseif(popNum==10)
				v_l=3;v_s=-3
			elseif(popNum==11)
				v_l=3;v_s=3
			elseif(popNum==12)
				v_l=3;v_s=2
			elseif(popNum==13)
				v_l=3;v_s=-2
			elseif(popNum==14)
				v_l=3;v_s=-1
			elseif(popNum==15)
				v_l=3;v_s=1
			elseif(popNum==16)
				v_l=3;v_s=-0
			endif
			break
		case -1: // control being killed
			break
	endswitch
	//print popNum
	return 0
End

Function Polaris_makewave(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	NVar v_pol_st=root:Polaris:misc:v_theta_s
	NVAR v_pol_del=root:Polaris:misc:v_polar_del
	NVAR v_pol_point=root:Polaris:misc:v_polar_point
	
	NVAR v_pol_en=root:Polaris:misc:v_theta_e
	
	NVAR v_tilt_st=root:Polaris:misc:v_tilt_s
	NVAR v_tilt_del=root:Polaris:misc:v_tilt_del
	NVAR v_tilt_point=root:Polaris:misc:v_tilt_point
	
	NVAR v_tilt_en=root:Polaris:misc:v_tilt_e
	
	NVAR v_Est=root:Polaris:misc:v_Est
	NVAR v_Edel=root:Polaris:misc:v_Edel
	NVAR v_Epoint=root:Polaris:misc:v_Epoint
	
	NVAR v_Een=root:Polaris:misc:v_Een

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			v_Een=v_Est+v_Edel*(v_Epoint-1)
			v_pol_en=v_pol_st+v_pol_del*(v_pol_point-1)
			v_tilt_en=v_tilt_st+v_tilt_del*(v_tilt_point-1)
			//print dval
			//ifchange end
			//v_Edel=(v_Een-v_Est)/(v_Epoint-1)
			//v_Edel=(v_Een-v_Est)/(v_Epoint-1)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Polaris_makewave_setend(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	NVar v_pol_st=root:Polaris:misc:v_theta_s
	NVAR v_pol_del=root:Polaris:misc:v_polar_del
	NVAR v_pol_point=root:Polaris:misc:v_polar_point
	
	NVAR v_pol_en=root:Polaris:misc:v_theta_e
	
	NVAR v_tilt_st=root:Polaris:misc:v_tilt_s
	NVAR v_tilt_del=root:Polaris:misc:v_tilt_del
	NVAR v_tilt_point=root:Polaris:misc:v_tilt_point
	
	NVAR v_tilt_en=root:Polaris:misc:v_tilt_e
	
	NVAR v_Est=root:Polaris:misc:v_Est
	NVAR v_Edel=root:Polaris:misc:v_Edel
	NVAR v_Epoint=root:Polaris:misc:v_Epoint
	
	NVAR v_Een=root:Polaris:misc:v_Een

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			//v_Een=v_Est+v_Edel*(v_Epoint-1)
			//v_pol_en=v_pol_st+v_pol_del*(v_pol_point-1)
			//v_tilt_en=v_tilt_st+v_tilt_del*(v_tilt_point-1)
			//print dval
			//ifchange end
			v_Edel=(v_Een-v_Est)/(v_Epoint-1)
			v_pol_del=(v_pol_en-v_pol_st)/(v_pol_point-1)
			v_tilt_del=(v_tilt_en-v_tilt_st)/(v_tilt_point-1)
			//v_Edel=(v_Een-v_Est)/(v_Epoint-1)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End




function Polaris_makewave_check()
	NVar v_pol_st=root:Polaris:misc:v_theta_s
	NVAR v_pol_del=root:Polaris:misc:v_polar_del
	NVAR v_pol_point=root:Polaris:misc:v_polar_point
	
	variable v_pol_st_pre, v_pol_del_pre, v_pol_point_pre
	
	//NVAR v_pol_en=root:Polaris:misc:v_theta_e
	
	NVAR v_tilt_st=root:Polaris:misc:v_tilt_s
	NVAR v_tilt_del=root:Polaris:misc:v_tilt_del
	NVAR v_tilt_point=root:Polaris:misc:v_tilt_point
	
	variable v_tilt_st_pre, v_tilt_del_pre, v_tilt_point_pre
	
	//NVAR v_tilt_en=root:Polaris:misc:v_tilt_e
	
	NVAR v_E_st=root:Polaris:misc:v_Est
	NVAR v_E_del=root:Polaris:misc:v_Edel
	NVAR v_E_point=root:Polaris:misc:v_Epoint
	
	variable 	v_E_st_pre, v_E_del_pre, v_E_point_pre
	
	NVAR v_kind=root:Polaris:misc:v_kind
	
	//NVAR v_Een=root:Polaris:misc:v_Een

	SetDataFolder root:Polaris:outuput
	wave/Z pre_w
	variable make_1or_0=0
	if(waveExists(pre_w)==1)
		if (v_kind==1)
			v_pol_st_pre=DimOffset(pre_w,1)
			v_pol_del_pre=DimDelta(pre_w,1)
			v_pol_point_pre=DimSize(pre_w,1)
					
			v_E_st_pre=DimOffset(pre_w,0)
			v_E_del_pre=DimDelta(pre_w,0)
			v_E_point_pre=DimSize(pre_w,0)
					
			if(v_pol_st==v_pol_st_pre && v_pol_del==v_pol_del_pre && v_pol_point==v_pol_point_pre\
			 && v_E_st==v_E_st_pre && v_E_del==v_E_del_pre && v_E_point==v_E_point_pre\
			  && wavedims(pre_w)==2)
				make_1or_0=0
			else
				make_1or_0=1
			endif
		elseif(v_kind==2)
			v_pol_st_pre=DimOffset(pre_w,0)
			v_pol_del_pre=DimDelta(pre_w,0)
			v_pol_point_pre=DimSize(pre_w,0)
			
			v_tilt_st_pre=DimOffset(pre_w,1)
			v_tilt_del_pre=DimDelta(pre_w,1)
			v_tilt_point_pre=DimSize(pre_w,1)
			if(v_pol_st==v_pol_st_pre && v_pol_del==v_pol_del_pre && v_pol_point==v_pol_point_pre\
			 && v_tilt_st==v_tilt_st_pre && v_tilt_del==v_tilt_del_pre && v_tilt_point==v_tilt_point_pre\
			 && wavedims(pre_w)==2)
				make_1or_0=0
			else
				make_1or_0=1
			endif
		elseif(v_kind==3)
			v_pol_st_pre=DimOffset(pre_w,0)
			v_pol_del_pre=DimDelta(pre_w,0)
			v_pol_point_pre=DimSize(pre_w,0)
			
			v_tilt_st_pre=DimOffset(pre_w,1)
			v_tilt_del_pre=DimDelta(pre_w,1)
			v_tilt_point_pre=DimSize(pre_w,1)
					
			v_E_st_pre=DimOffset(pre_w,0)
			v_E_del_pre=DimDelta(pre_w,0)
			v_E_point_pre=DimSize(pre_w,0)
			
			if(v_pol_st==v_pol_st_pre && v_pol_del==v_pol_del_pre && v_pol_point==v_pol_point_pre\
			 && v_tilt_st==v_tilt_st_pre && v_tilt_del==v_tilt_del_pre && v_tilt_point==v_tilt_point_pre\
			 && v_E_st==v_E_st_pre && v_E_del==v_E_del_pre && v_E_point==v_E_point_pre\
			  && wavedims(pre_w)==3)
			   	make_1or_0=0
			else
				make_1or_0=1
			endif
		endif
				
	else
		make_1or_0=1
	endif
			
	if(make_1or_0==1)
		PrepareWave(v_kind, v_E_point,v_E_st,v_E_del, v_pol_point, v_pol_st, v_pol_del,v_tilt_point, v_tilt_st, v_tilt_del)
	endif
End




Function Polaris_tiltmap2DA30_SliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	NVAR v_tiltamp2DA30=root:Polaris:misc:tiltmap2DA30
	NVAR v_check_angle_or_k=root:Polaris:misc:v_check_angle_or_k

	switch( sa.eventCode )
		case -3: // Control received keyboard focus
		case -2: // Control lost keyboard focus
		case -1: // Control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				if(curval==0)
					v_tiltamp2DA30=1 ;v_check_angle_or_k=0
				elseif(curval==1)
					v_tiltamp2DA30=1 ;v_check_angle_or_k=1
				elseif(curval==2)
					v_tiltamp2DA30=0 ;v_check_angle_or_k=0
				elseif(curval==3)
					v_tiltamp2DA30=0 ;v_check_angle_or_k=1
				endif
			endif
			break
	endswitch
	 //print curval
	return 0
End




Function Polaris_AutoDisplayStyle(w)
	wave w
	NVAR v_kind=root:Polaris:misc:v_kind
	NVAR v_check_angle_or_k=root:Polaris:misc:v_check_angle_or_k
	if(v_kind==2 && v_check_angle_or_k==1)
		ModifyGraph/Z swapXY=0
		Label/Z left "\\f02k\\By(⊥)\\M\\f00 (Å\\S-1\\M)"
		Label/Z bottom "\\f02k\\Bx(∥)\\M\\f00 (Å\\S-1\\M)"
	elseif(v_kind==2 && v_check_angle_or_k==0)
		ModifyGraph/Z swapXY=0
		Label/Z left "\\f02ψ\f00 (deg.)"
		Label/Z bottom "\\f02θ\f00 (deg.)"
	elseif(v_kind==1 && v_check_angle_or_k==0)
		ModifyGraph/Z swapXY=1
		Label/Z left "\\f02E - E\\B\\f00F\\M (eV)"
		Label/Z bottom "\\f02θ\f00 (deg.)"
	elseif(v_kind==1 && v_check_angle_or_k==1)
		ModifyGraph/Z swapXY=1
		Label/Z left "\\f02E - E\\B\\f00F\\M (eV)"
		Label/Z bottom "\\f02k\\Bx(∥)\\M\\f00 (Å\\S-1\\M)"
		
	endif
		ModifyGraph/Z margin(right)=2,margin(top)=2, tlOffset(left)=-2,tlOffset(bottom)=-4, zero(left)=0
	
	ModifyGraph/Z tick=2, ZisZ=0, standoff=0, tickUnit=1, margin(left)=21, margin(bottom)=21
	GetAxis/Q right
	if (V_flag)
		ModifyGraph/Z mirror(left)=1
	endif
	GetAxis/Q top
	if (V_flag)
		ModifyGraph/Z mirror(bottom)=1
	endif
	AutoPositionWindow/E
End


Function Polaris_Button_plot(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	NVAR v_kind=root:Polaris:misc:v_kind
	NVAR v_plot=root:Polaris:misc:v_select_display
	string pol,tot,select
	pol="polwav"
	SVAR cubspfunc=root:Polaris:misc:st_orb
	tot=cubspfunc+"_tot_mat"
	if(v_plot==1)//pol
		wave cubsp=$pol
		select=pol
	elseif(v_plot==2)//orb
		wave cubsp=$cubspfunc	
		select=cubspfunc
	elseif(v_plot==3)//all
		wave cubsp=$tot	
		select=tot
	endif
	

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			///plot
			SetDataFolder root:Polaris:outuput
			if(WinType(select+"0")==0)
					variable dim,ni=0
					dim=WaveDims(cubsp)
					switch(dim)			
						case 3:	//volume
							if(exists("CrossSectionViewer")==6)
								//CrossSectionViewer(cubsp)
							endif
						break
						case 2:	//image
						case 1:	// trace
							if (v_kind==2)
								Display/N=$select/K=1/W=(44,44,250,250)
								//Dowindow/C $select
							elseif(v_kind==1)
								Display/N=$select/K=1/W=(1,44,198,304)
							endif
							if (dim==2)	//image
								AppendImage cubsp
								ModifyImage ''#(ni) ctab= {,,BlueHot256,0}, ctabAutoscale=1, interpolate=-1
								ni+=1
							endif
					endswitch
					///if (!new)
						Polaris_AutoDisplayStyle(cubsp)
					//endif
				break
			endif
			Polaris_AutoDisplayStyle(cubsp)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Polaris_ButtonProc_output(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			Button Output title="\\f03Busy"
			Polaris_makewave_check()
			Polaris_cal()
			Button Output title="\\f03Output"
			
			///plot
			
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Polaris_paraset(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	NVAR v_l=root:Polaris:misc:v_l
	NVAR v_m=root:Polaris:misc:v_s
	NVAR v_kind=root:Polaris:misc:v_kind
	string tot
	SVAR cubspfunc=root:Polaris:misc:st_orb
	tot=cubspfunc+"_tot_mat"
	wave/Z cubsp=$cubspfunc
	wave/Z polwav
	wave/Z Total_mat_el=$tot
	
	//励起光,仕事関数, 内部ポテンシャル
	NVAR v_hn=root:Polaris:misc:v_hn
	NVAR v_W=root:Polaris:misc:v_W
	NVAR v_V0=root:Polaris:misc:v_V0
	NVAR v_escdep=root:Polaris:misc:v_escdep
	NVAR v_escdep_fit=root:Polaris:misc:v_escdep_fit

	//Wave種類
	NVAR v_kind=root:Polaris:misc:v_kind
	//アジマス角
	NVAR v_azi=root:Polaris:misc:v_azi
	NVAR v_tilt=root:Polaris:misc:v_tilt
	
	//lightpol
	NVAR v_Lpol=root:Polaris:misc:light_pol:v_guzal
	NVAR v_Cpol=root:Polaris:misc:light_pol:v_Delta
	NVAR v_incdeg_pol=root:Polaris:misc:light_pol:v_inc_angle
	NVAR v_Intofs=root:Polaris:misc:v_Intofs
	
	NVAR v_liveupdate=root:Polaris:misc:v_liveupdate
	NVAR v_tiltamp2DA30=root:Polaris:misc:tiltmap2DA30
	NVAR v_check_angle_or_k=root:Polaris:misc:v_check_angle_or_k
	

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			v_escdep_fit=(143/(v_hn-v_w)^2+0.054*(v_hn-v_w)^0.5)*10
			if(v_liveupdate==1)
				Button Output title="\\f03Busy"
				Polaris_cal()
				Button Output title="\\f03Output"
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Polaris_nan2zero(x)
	variable/D x
	return (numtype(x) ? 0 : x)
End

Function POLARIS_MACRO() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /N=POLARIS_window /K=1 /W=(329,277,1009,896) as "POLARIS_panel"
	ModifyPanel cbRGB=(49151,65535,49151), fixedSize=1
	SetDrawLayer UserBack
	SetDrawEnv linethick= 0,fillfgc= (51664,44236,58982)
	DrawRect -17,194,739,363
	SetDrawEnv linefgc= (1,65535,33232),fillfgc= (1,65535,33232)
	SetDrawEnv save
	SetDrawEnv linefgc= (0,0,0),arrow= 2
	DrawLine 410,42,410,157
	SetDrawEnv linefgc= (0,0,0),arrow= 1
	DrawLine 354,92,465,109
	SetDrawEnv fsize= 20
	DrawText 475,134,"\\f02x"
	SetDrawEnv fsize= 20
	DrawText 435,51,"\\f02y"
	SetDrawEnv linefgc= (0,0,0),arrow= 1
	DrawLine 461,79,296,151
	SetDrawEnv fsize= 16
	DrawText 32,26,"POlarization dependent ARpes and Incident light dependent  Simulation (Polaris) Macro"
	SetDrawEnv fsize= 15
	DrawText 312,395,"Light polarization Term"
	DrawText 101,236,"Start"
	SetDrawEnv fsize= 15
	DrawText 32,51,"Orbital selection Term"
	SetDrawEnv fsize= 15
	DrawText 32,100,"Kinds of Analyzer Term"
	DrawText 62,190,"DA30 (angle)"
	DrawText 62,146,"Tilt map (angle)"
	SetDrawEnv fsize= 15
	DrawText 433,217,"Kinds of wave dimenssion "
	SetDrawEnv fsize= 15
	DrawText 32,395,"Chage parameter"
	SetDrawEnv fsize= 15
	DrawText 184,176,"Light pol. def."
	DrawText 268,647,"Output data"
	SetDrawEnv fsize= 15
	DrawText 216,195,"\\f03ε = \\f00\\f02ε \\f00( e\\S\\f02iδ\\M\\f00cos\\f02α\\f00 sin\\f02ξ\\f00, cos\\f02ξ\\f00cos\\f02ψ\\f00, e\\S\\f02iδ\\M\\f00sin\\f02α\\f00 sin\\f02ξ\\f00 )\\St"
	DrawText 100,100,"\\f00"
	SetDrawEnv fsize= 20
	DrawText 310,166,"\\f02z"
	DrawText 333,263,"deg. / Å\\S-1"
	DrawText 334,289,"deg. / Å\\S-1"
	DrawText 146,503,"deg. "
	DrawText 62,125,"Tilt map (Momentum)"
	DrawText 62,169,"DA30 (Momentum)"
	DrawText 446,283,"\\f02E - k\\Bx\\M\\f00 or \\f02E - θ"
	DrawText 446,261,"\\f02k\\Bx\\M - k\\By\\M\\f00 or \\f02θ- ψ"
	DrawText 446,241,"\\f02E (k\\Bx\\M, k\\By\\M\\f00) or \\f02E (θ, ψ)"
	DrawText 159,236,"Step"
	DrawText 201,236,"Data point"
	SetDrawEnv linethick= 0.8,linefgc= (1,4,52428),linecap= 2,arrow= 3
	DrawLine 265,152,265,99
	SetDrawEnv linefgc= (65535,0,0)
	DrawLine 453,107,387,134
	SetDrawEnv linefgc= (65535,0,0)
	DrawLine 375,93,309,120
	SetDrawEnv linefgc= (65535,0,0)
	DrawLine 389,134,310,121
	SetDrawEnv linethick= 0.8,linefgc= (1,4,52428),linecap= 1,arrow= 1
	DrawLine 218,139,409,100
	SetDrawEnv fsize= 20
	DrawText 247,109,"\\f03ε"
	SetDrawEnv linefgc= (0,0,0),dash= 1,fillpat= 0
	DrawOval 253,99,276,154
	SetDrawEnv linefgc= (0,0,0),linecap= 2,arrow= 3,arrowlen= 5,arrowfat= 1
	DrawLine 253,131,276,126
	SetDrawEnv linefgc= (0,0,0),arrow= 2,fillpat= 0
	DrawArc  264,118,26.0192236625154,0,91.0230301886679
	SetDrawEnv fsize= 20
	DrawText 280,100,"\\f02ξ"
	SetDrawEnv linefgc= (0,0,0),arrow= 2,arrowlen= 5,fillpat= 0
	DrawArc  381,90,26.0192236625154,-135.881403996582,-101.496563017586
	SetDrawEnv linefgc= (0,0,0),dash= 3,arrow= 1
	DrawLine 317,94,365,115
	SetDrawEnv fsize= 20
	DrawText 305,102,"\\f02α"
	SetDrawEnv linefgc= (0,0,0),arrow= 2,arrowlen= 5,fillpat= 0
	DrawArc  462,110,15,89.2152,-155.323
	SetDrawEnv linefgc= (0,0,0),arrow= 1,arrowlen= 5,fillpat= 0
	DrawArc  411,51,15,-168.870810710389,-39.4400527366905
	DrawText 465,102,"\\f02ψ\\f00: Tilt angle"
	DrawText 426,70,"\\f02θ\\f00: Polar angle"
	SetDrawEnv linefgc= (0,0,0),arrow= 1,arrowlen= 5,fillpat= 0
	DrawArc  317,144,13,-15.5241,105.255118703058
	DrawText 330,165,"\\f02φ\\f00: Azimuth angle"
	SetDrawEnv fsize= 15
	DrawText 32,216,"Make wave Term"
	DrawText 284,236,"End"
	DrawText 497,580,"Note:\rξ = 0 and δ = 0 => s-pol\rξ = 90 and δ = 0 => p-pol\rξ = 45 and δ = 90 => right cir.-pol\rξ = 45 and δ = -90 => left cir.-pol"
	SetDrawEnv fsize= 15
	DrawText 32,483,"Orbital Term"
	DrawText 57,526,"(only \\f02E- k\\Bx\\M\\f00 or \\f02E- θ\\f00 )"
	DrawText 148,464,"=1430/\\f02E\\S2\\M\\Bkin\\M\\f00+0.54\\f02√E\\Bkin\\M\\f00(Cf.)"
	DrawText 556,485,"only pol. term"
	DrawText 556,464,"only orb. term"
	DrawText 556,442,"pol. × orb. term"
	DrawText 148,546,"deg. / Å\\S-1\\M"
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 5,267,325,289
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0),fillbgc= (44253,29492,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0,0)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillfgc= (0,65535,0)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= 0,fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linethick= 0,linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linethick= 0,linefgc= (51664,44236,58982),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 1,267,329,289
	SetDrawEnv linethick= 0,linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,242,326,264
	SetDrawEnv linefgc= (0,65535,0),fillpat= -1,fillbgc= (51664,44236,58982)
	DrawRect 6,267,326,289
	SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
	DrawRect 6,294,326,316
	SetVariable v_pol_st,pos={16.00,244.00},size={122.00,20.00},bodyWidth=40,proc=Polaris_makewave,title="Polar\\f02k\\Bx\\f00\\M range"
	SetVariable v_pol_st,fSize=12
	SetVariable v_pol_st,limits={-inf,inf,0},value= root:Polaris:misc:v_theta_s,live= 1
	SetVariable v_tilt_st,pos={23.00,267.00},size={115.00,20.00},bodyWidth=40,proc=Polaris_makewave,title="Tilt/\\f02k\\By\\f00\\M range"
	SetVariable v_tilt_st,fSize=12
	SetVariable v_tilt_st,limits={-inf,inf,0},value= root:Polaris:misc:v_tilt_s,live= 1
	SetVariable v_azi,pos={57.00,484.00},size={90.00,18.00},bodyWidth=40,proc=Polaris_paraset,title="Azimuth"
	SetVariable v_azi,fSize=12
	SetVariable v_azi,limits={-inf,inf,0},value= root:Polaris:misc:v_azi,live= 1
	Slider ang2wav,pos={39.00,109.00},size={21.00,82.00},proc=Polaris_tiltmap2DA30_SliderProc
	Slider ang2wav,limits={0,3,1},value= 3,ticks= 0
	PopupMenu angular_momentum_l,pos={53.00,50.00},size={100.00,23.00},bodyWidth=100,proc=Polaris_PopMenu_orbital_selection_Proc
	PopupMenu angular_momentum_l,fSize=12
	PopupMenu angular_momentum_l,mode=2,popvalue="px",value= #"\"s;px;py;pz;dxy;dx^2-y^2;dxz;dyz;dz^2;fx(x^2-3y^2);fy(3x^2-y^2);fxyz;fz(x^2-y^2);fxz^2;fyz^2;fz^3\""
	Slider ang2wav1,pos={427.00,222.00},size={21.00,60.00},proc=Polaris_SliderProc_set_wavedim
	Slider ang2wav1,limits={1,3,1},variable= root:Polaris:misc:v_kind,ticks= 0
	SetVariable v_pol_end1,pos={153.00,269.00},size={40.00,18.00},bodyWidth=40,proc=Pol_ARPES_Set_V_tilt_end
	SetVariable v_pol_end1,fSize=12,limits={-inf,inf,0}
	SetVariable v_en_st,pos={35.00,295.00},size={103.00,20.00},bodyWidth=40,proc=Polaris_makewave,title="\\f02E - E\\f00\\BF\\M (eV)"
	SetVariable v_en_st,fSize=12
	SetVariable v_en_st,limits={-inf,inf,0},value= root:Polaris:misc:v_Est,live= 1
	SetVariable v_en_del,pos={155.00,295.00},size={40.00,18.00},bodyWidth=40,proc=Polaris_makewave,title=" "
	SetVariable v_en_del,fSize=12
	SetVariable v_en_del,limits={-inf,inf,0},value= root:Polaris:misc:v_Edel,live= 1
	SetVariable v_polardel,pos={155.00,244.00},size={40.00,18.00},bodyWidth=40,proc=Polaris_makewave,title=" "
	SetVariable v_polardel,fSize=12
	SetVariable v_polardel,limits={-inf,inf,0},value= root:Polaris:misc:v_polar_del,live= 1
	SetVariable v_tilt_del,pos={155.00,269.00},size={40.00,18.00},bodyWidth=40,proc=Polaris_makewave,title=" "
	SetVariable v_tilt_del,fSize=12
	SetVariable v_tilt_del,limits={-inf,inf,0},value= root:Polaris:misc:v_tilt_del,live= 1
	SetVariable v_en_point,pos={203.00,295.00},size={50.00,18.00},bodyWidth=50,proc=Polaris_makewave,title=" "
	SetVariable v_en_point,fSize=12
	SetVariable v_en_point,limits={1,inf,1},value= root:Polaris:misc:v_Epoint,live= 1
	SetVariable v_pol_point,pos={203.00,244.00},size={50.00,18.00},bodyWidth=50,proc=Polaris_makewave,title=" "
	SetVariable v_pol_point,fSize=12
	SetVariable v_pol_point,limits={1,inf,1},value= root:Polaris:misc:v_polar_point,live= 1
	SetVariable v_tilt_point,pos={203.00,269.00},size={50.00,18.00},bodyWidth=50,proc=Polaris_makewave,title=" "
	SetVariable v_tilt_point,fSize=12
	SetVariable v_tilt_point,limits={1,inf,1},value= root:Polaris:misc:v_tilt_point,live= 1
	SetVariable v_Lpol,pos={323.00,409.00},size={108.00,18.00},bodyWidth=60,proc=Polaris_paraset,title="\\f02ξ\\f00 (deg.)"
	SetVariable v_Lpol,fSize=12
	SetVariable v_Lpol,limits={-inf,inf,0},value= root:Polaris:misc:light_pol:v_guzal
	SetVariable v_Cpol,pos={323.00,436.00},size={151.00,18.00},bodyWidth=60,proc=Polaris_paraset,title="\\f02δ\\f00 (deg.): Phase"
	SetVariable v_Cpol,fSize=12
	SetVariable v_Cpol,limits={-inf,inf,0},value= root:Polaris:misc:light_pol:v_Delta
	SetVariable v_Lpol2,pos={323.00,457.00},size={196.00,18.00},bodyWidth=60,proc=Polaris_paraset,title="\\f02α\\f00 (deg.): incident angle"
	SetVariable v_Lpol2,fSize=12
	SetVariable v_Lpol2,limits={-90,90,0},value= root:Polaris:misc:light_pol:v_inc_angle,live= 1
	SetVariable v_hn,pos={57.00,400.00},size={85.00,18.00},bodyWidth=40,proc=Polaris_paraset,title="\\f02hn\\f00 (eV)"
	SetVariable v_hn,fSize=12,limits={0,inf,0},value= root:Polaris:misc:v_hn,live= 1
	SetVariable v_W0,pos={56.00,422.00},size={168.00,18.00},bodyWidth=40,proc=Polaris_paraset,title="W (eV): Work function"
	SetVariable v_W0,fSize=12,limits={0,inf,0},value= root:Polaris:misc:v_W,live= 1
	SetVariable v_pol_en,pos={271.00,244.00},size={50.00,18.00},proc=Polaris_makewave_setend,title=" "
	SetVariable v_pol_en,fSize=12
	SetVariable v_pol_en,limits={-inf,inf,0},value= root:Polaris:misc:v_theta_e
	SetVariable v_tilt_en,pos={271.00,270.00},size={50.00,18.00},proc=Polaris_makewave_setend,title=" "
	SetVariable v_tilt_en,fSize=12
	SetVariable v_tilt_en,limits={-inf,inf,0},value= root:Polaris:misc:v_tilt_e
	SetVariable v_E_en,pos={271.00,295.00},size={50.00,18.00},proc=Polaris_makewave_setend,title=" "
	SetVariable v_E_en,fSize=12,limits={-inf,inf,0},value= root:Polaris:misc:v_Een
	Button output,pos={544.00,373.00},size={100.00,20.00},proc=Polaris_ButtonProc_output,title="\\f03Output"
	CheckBox liveupdate,pos={164.00,377.00},size={86.00,16.00},title=" Live update"
	CheckBox liveupdate,fSize=12,variable= root:Polaris:misc:v_liveupdate
	SetVariable v_tilt,pos={57.00,527.00},size={82.00,20.00},bodyWidth=40,proc=Polaris_paraset,title="Tilt/\\f02k\\By\\f00\\M "
	SetVariable v_tilt,fSize=12
	SetVariable v_tilt,limits={-90,90,0},value= root:Polaris:misc:v_tilt,live= 1
	Button plot,pos={584.00,490.00},size={50.00,20.00},proc=Polaris_Button_plot,title="\\f03plot"
	SetVariable v_V0,pos={323.00,482.00},size={172.00,20.00},bodyWidth=40,proc=Polaris_paraset,title="V\\B0\\M (eV): Inner potential"
	SetVariable v_V0,fSize=12
	SetVariable v_V0,limits={-inf,inf,0},value= root:Polaris:misc:v_V0,live= 1
	SetVariable v_escdep,pos={323.00,509.00},size={157.00,18.00},bodyWidth=40,proc=Polaris_paraset,title="λ (Å): Probing depth"
	SetVariable v_escdep,fSize=12
	SetVariable v_escdep,limits={1e-05,inf,0},value= root:Polaris:misc:v_escdep,live= 1
	SetVariable v_offset,pos={323.00,535.00},size={129.00,18.00},bodyWidth=40,proc=Polaris_paraset,title="Intensity offset"
	SetVariable v_offset,fSize=12
	SetVariable v_offset,limits={-inf,inf,0},value= root:Polaris:misc:v_Intofs,live= 1
	Slider select_plot,pos={530.00,425.00},size={21.00,60.00}
	Slider select_plot,limits={1,3,1},variable= root:Polaris:misc:v_select_display,ticks= 0
	SetVariable v_probing_dep,pos={57.00,443.00},size={88.00,20.00},bodyWidth=50,proc=Polaris_paraset,title="λ\\Bfit\\M (Å)"
	SetVariable v_probing_dep,fSize=12,frame=0
	SetVariable v_probing_dep,limits={-inf,inf,0},value= root:Polaris:misc:v_escdep_fit,live= 1
	TitleBox xslice,pos={518.00,3.00},size={50.00,20.00}
	Button button_setscaleformwave,pos={163.00,202.00},size={100.00,20.00},proc=Polaris_ButtonProc_setscalefromwave,title="Set from Wave"
	//SetWindow kwTopWin,hook(CSwin)=CSwinhook
EndMacro

Function Polaris_SliderProc_set_wavedim(sa) : SliderControl
	STRUCT WMSliderAction &sa
	NVAR v_kind=root:Polaris:misc:v_kind

	switch( sa.eventCode )
		case -3: // Control received keyboard focus
		case -2: // Control lost keyboard focus
		case -1: // Control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				if(v_kind==1)
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					DrawRect 6,242,326,264
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					SetDrawEnv fillpat= -1
					DrawRect 6,267,326,289
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					DrawRect 6,294,326,316
					v_kind=1
				elseif(v_kind==2)
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					DrawRect 6,242,326,264
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					DrawRect 6,267,326,289
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					SetDrawEnv fillpat= -1
					DrawRect 6,294,326,316
					v_kind=2
				elseif(v_kind==3)
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					DrawRect 6,242,326,264
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					DrawRect 6,267,326,289
					SetDrawEnv linefgc= (0,65535,0),fillbgc= (51664,44236,58982)
					DrawRect 6,294,326,316
					v_kind=3
				endif
			endif
			break
	endswitch

	return 0
End

Function Polaris_set_scale_from_wave()
	string tws=GetBrowserSelection(0)
	NVar v_pol_st=root:Polaris:misc:v_theta_s // start
	NVAR v_pol_del=root:Polaris:misc:v_polar_del // step
	NVAR v_pol_point=root:Polaris:misc:v_polar_point // points
	
	NVAR v_theta_e=root:Polaris:misc:v_theta_e
	
	NVAR v_tilt_st=root:Polaris:misc:v_tilt_s
	NVAR v_tilt_del=root:Polaris:misc:v_tilt_del
	NVAR v_tilt_point=root:Polaris:misc:v_tilt_point
	
	NVAR v_tilt_e=root:Polaris:misc:v_tilt_e
	
	NVAR v_Est=root:Polaris:misc:v_Est
	NVAR v_Edel=root:Polaris:misc:v_Edel
	NVAR v_Epoint=root:Polaris:misc:v_Epoint
	
	NVAR v_Een=root:Polaris:misc:v_Een
	
	NVAR v_kind=root:Polaris:misc:v_kind // 1:E-po, 2:po-ti, 3:E-po-ti
	
	variable go=0,preset=0
	
	if(WaveExists($tws))
		if(WaveDims($tws)==2)
			if(v_kind==1 || v_kind==2)
				go=1
			else
				go=0
				print "Select 3D Wave"
			endif
		elseif(WaveDims($tws)==3)
			if(v_kind==3)
				go=1
			else
				go=0
				print "Select 2D Wave"
			endif
		else
		 go=0
		endif
	else
		go=0
		preset=1
		print "Select Wave"
	endif
	
	if(go)
		wave tw=$tws
		if(v_kind==1)// 1:E-po, 2:po-ti, 3:E-po-ti
			v_Est=DimOffset(tw,0)
			v_Edel=DimDelta(tw,0)
			v_Epoint=DimSize(tw,0)
			
			v_Een=v_Est+v_Edel*(v_Epoint-1)
			
			v_pol_st=DimOffset(tw,1)
			v_pol_del=DimDelta(tw,1)
			v_pol_point=DimSize(tw,1)
			
			v_theta_e=v_pol_st+v_pol_del*(v_pol_point-1)
		elseif(v_kind==2)
			v_pol_st=DimOffset(tw,0)
			v_pol_del=DimDelta(tw,0)
			v_pol_point=DimSize(tw,0)
			
			v_theta_e=v_pol_st+v_pol_del*(v_pol_point-1)
			
			v_tilt_st=DimOffset(tw,1)
			v_tilt_del=DimDelta(tw,1)
			v_tilt_point=DimSize(tw,1)
			
			v_tilt_e=v_tilt_st+v_tilt_del*(v_tilt_point-1)
		elseif(v_kind==3)
			v_Est=DimOffset(tw,0)
			v_Edel=DimDelta(tw,0)
			v_Epoint=DimSize(tw,0)
			
			v_Een=v_Est+v_Edel*(v_Epoint-1)
			
			v_pol_st=DimOffset(tw,1)
			v_pol_del=DimDelta(tw,1)
			v_pol_point=DimSize(tw,1)
			
			v_theta_e=v_pol_st+v_pol_del*(v_pol_point-1)
			
			v_tilt_st=DimOffset(tw,2)
			v_tilt_del=DimDelta(tw,2)
			v_tilt_point=DimSize(tw,2)
			
			v_tilt_e=v_tilt_st+v_tilt_del*(v_tilt_point-1)
		endif
	endif
	
	if(preset==1)
		v_Est=-2
		v_Edel=0.01
		v_Epoint=301
		
		v_Een=v_Est+v_Edel*(v_Epoint-1)
		
		v_pol_st=-15
		v_pol_del=0.1
		v_pol_point=301
		
		v_theta_e=v_pol_st+v_pol_del*(v_pol_point-1)
		
		v_tilt_st=-15
		v_tilt_del=0.1
		v_tilt_point=301
		
		v_tilt_e=v_tilt_st+v_tilt_del*(v_tilt_point-1)	
	
	endif
End

Function Polaris_ButtonProc_setscalefromwave(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			Polaris_set_scale_from_wave()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function Polaris_cal()
	SetDataFolder root:Polaris:outuput
				//代入兼プロット
	NVAR v_tiltamp2DA30=root:Polaris:misc:tiltmap2DA30
	NVAR v_check_angle_or_k=root:Polaris:misc:v_check_angle_or_k
	
	NVAR v_l=root:Polaris:misc:v_l
	NVAR v_m=root:Polaris:misc:v_s
	//励起光,仕事関数, 内部ポテンシャル
	NVAR v_hn=root:Polaris:misc:v_hn
	NVAR v_W=root:Polaris:misc:v_W
	NVAR v_V0=root:Polaris:misc:v_V0
	NVAR v_escdep=root:Polaris:misc:v_escdep
	//Wave種類
	NVAR v_kind=root:Polaris:misc:v_kind
	//アジマス角
	NVAR v_azi=root:Polaris:misc:v_azi
	NVAR v_tilt=root:Polaris:misc:v_tilt
	
	//lightpol
	NVAR v_Lpol=root:Polaris:misc:light_pol:v_guzal
	NVAR v_Cpol=root:Polaris:misc:light_pol:v_Delta
	NVAR v_incdeg_pol=root:Polaris:misc:light_pol:v_inc_angle
	NVAR v_Intofs=root:Polaris:misc:v_Intofs
			
	string l_s,m_s,tot
	SVAR cubspfunc=root:Polaris:misc:st_orb
	wave pre_w
	if(v_l==0)
		l_s="s"
		m_s=""
	elseif(v_l==1)
		l_s="p"
		if(v_m==1)
			m_s="x"
		elseif(v_m==0)
			m_s="z"
		elseif(v_m==-1)
			m_s="y"
		endif
	elseif(v_l==2)
		l_s="d"
		if(v_m==2)
			m_s="x2_y2"
		elseif(v_m==1)
			m_s="xz"
		elseif(v_m==0)
			m_s="z2"
		elseif(v_m==-1)
			m_s="yz"
		elseif(v_m==-2)
			m_s="xy"
		endif
	elseif(v_l==3)
		l_s="f"
		if(v_m==3)
			m_s="y(3x2_y2)"	
		elseif(v_m==2)
			m_s="xyz"
		elseif(v_m==1)
			m_s="yz2"
		elseif(v_m==0)
			m_s="z3"
		elseif(v_m==-1)
			m_s="xz2"
		elseif(v_m==-2)
			m_s="z(x2_y2)"
		elseif(v_m==-3)
			m_s="x(x2_3y2)"
		endif
	endif
	
	cubspfunc=l_s+m_s+"orb"
	tot=cubspfunc+"_tot_mat"
	duplicate/o pre_w, polwav
	duplicate/o pre_w, $cubspfunc
	duplicate/o pre_w, $tot
	wave cubsp=$cubspfunc
	wave Total_mat_el=$tot
	//variable x,y,z
	if(v_tiltamp2DA30==0)//tiltmap
		if(v_kind==1)
			cubsp=(Polaris_Wigner_Rot_Matrix_tiltmap(v_l,v_m,v_hn,v_W, 0,v_tilt,v_azi,v_check_angle_or_k, v_kind, x,y,z))^2
			polwav=(Polaris_pol_light_tiltmap(v_Lpol, v_Cpol, v_incdeg_pol,v_V0, v_escdep, v_hn,v_W,x,y,v_tilt,v_kind, v_check_angle_or_k))^2+v_Intofs
			Total_mat_el=cubsp*polwav
		else
			cubsp=(Polaris_Wigner_Rot_Matrix_tiltmap(v_l,v_m,v_hn,v_W, 0,0,v_azi,v_check_angle_or_k, v_kind, x,y,z))^2
			polwav=(Polaris_pol_light_tiltmap(v_Lpol, v_Cpol, v_incdeg_pol,v_V0, v_escdep, v_hn,v_W,x,y,z,v_kind, v_check_angle_or_k))^2+v_Intofs
			Total_mat_el=cubsp*polwav
		 endif
				
			
	elseif(v_tiltamp2DA30==1)//DA30
		if(v_kind==1)
			cubsp=(Polaris_Wigner_Rot_Matrix_DA30(v_l,v_m,v_hn,v_W, 0,0,v_azi,v_check_angle_or_k, v_kind, x,y,v_tilt))^2
			polwav=(Polaris_pol_light_DA30(v_Lpol, v_Cpol, v_incdeg_pol,0,v_V0, v_escdep, v_hn,v_W,x,y,v_tilt,v_kind, v_check_angle_or_k))^2+v_Intofs
			Total_mat_el=cubsp*polwav
		else
			cubsp=(Polaris_Wigner_Rot_Matrix_DA30(v_l,v_m,v_hn,v_W, 0,0,v_azi,v_check_angle_or_k, v_kind, x,y,z))^2
			polwav=(Polaris_pol_light_DA30(v_Lpol, v_Cpol, v_incdeg_pol,0,v_V0, v_escdep, v_hn,v_W,x,y,z,v_kind, v_check_angle_or_k))^2+v_Intofs
			Total_mat_el=cubsp*polwav
		endif
	endif
end

