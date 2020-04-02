#include<stdio.h>
#include<math.h>
#define PI 3.14
double v,r,c,k,temp,h;

void Prop_oil(double temp);
void Prop_air(double temp);
void Prop_water(double temp);
double CALC_TCO(double e, double cc, double ch, double th, double tc);
double CALC_THO(double e, double cc, double ch, double th, double tc);
double MATERIAL_THERMALCONDUCTIVITY(int c);
double MATERIAL_VELOCITY(int a);
double TEMA_TUBESIDE_INNERDIA(int b);
double TEMA_TUBESIDE_OUTERDIA(int i);
double TEMA_ANNULARSIDE_INNERDIA(int a);
double TEMA_ANNULARSIDE_OUTERDIA(int c);
double CALC_RE(double r, double v, double d, double vi);
void makefile();
double CALC_PR(double vi, double cp, double k);
double CALC_NU(double rea, double pra);
double CALC_H(double n, double k, double d);
FILE *a1,*a2,*a3,*a4,*a5,*fp,*fp1,*fp2,*fp3,*fp4,*a6,*fp5,*fp6,*fp7,*fp8,*fp9;
double CALC_UC(double d0, double di, double hi, double k, double h0);
double CALC_UF(double u, double rf);
double CALC_CF(double uf, double uc);
double CALC_OD(double c);
double CALC_LMTD(double Tin, double Tout);
double calc_f(double re, double K, double D);
double calc_90(double a);
double calc_pressuredrop_Annulus(double f, double n, double l, double d, double vi, double rho, double h);
double calc_pressuredrop_Tubein(double f, double n, double l, double d, double u, double rho);

main()
{
	int x,y,z,z1,j,i,k1,N1;
	double D1,Pri,Pro,Tbh,Tbc,Tho,Tco,vmi,Nui,Nuo,Dei,Deo,hi,ho,Uc,Uf,CF,OD,Rft,Q,LMTD,Ao,Tcommin,Tcommout,L,N,fi,f0,ki,ko,P,Pin,Pout,N2,b1,b2,n,o;
	double Mh,Mc,Thi,Tci,Epsilon,rho_h,rho_c,Cph,Cpc,Ch,Cc,Cmin,Cmax,R,vis_h,vis_c,k_h,k_c,k_m,d0,di,ti,t0,Di,D0,Dh,De,Vmi,Vmo,Vmax,Rei,Reo,Aci,Aco;
	makefile();
	fscanf(fp,"%*s%lf %*s %*s%lf %*s %*s%lf %*s %*s%lf %*s %*s%lf %*s%lf %*s",&Mh,&Mc,&Thi,&Tci,&Epsilon,&L);
	fprintf(fp1,"Mass flow rate of hot fluid = %lf Kg/s. \n",Mh);
	fprintf(fp2,"Mass flow rate of cold fluid = %lf Kg/s. \n",Mc);
	fprintf(fp1,"Hot Fluid initial Temperature = %lf 'c . \n",Thi);
	fprintf(fp2,"Cold Fluid initial Temperature = %lf 'c . \n",Tci);
	fclose(fp);
	ENTER_FLUID_HOT:
	printf("Choose the Fluid flow as Hot Fluid : \n");
	printf("1. Water \n");
	printf("2. Air \n");
	printf("3. oil \n");
	printf("4. Other \n");
	printf("Enter the selected fluid : ");
	scanf("%d",&x);
	if(x > 4 || x < 1)
	{
		printf("Invalid option Selected. \n");
		goto ENTER_FLUID_HOT;
	}
	if(x == 1)
	{
		Prop_water(Thi);
		Cph = c;
	}
	if(x == 2)
	{
		Prop_air(Thi);
		Cph = c;
	}
	if(x == 3)
	{
		Prop_oil(Thi);
		Cph = c;
	}
	if(x == 4)
	{
		printf("Enter the Specific Heat of Hot Fluid at Temp %lf 'c (in J/Kg.K): ",Thi);
		scanf("%lf",&Cph);
	}
	ENTER_FLUID_COLD:
	printf("Choose the Fluid flow as Cold Fluid : \n");
	printf("1. Water \n");
	printf("2. Air \n");
	printf("3. oil \n");
	printf("4. Other \n");
	printf("Enter the selected fluid : ");
	scanf("%d",&y);
	if(y > 4 || y < 1)
	{
		printf("Invalid option Selected. \n");
		goto ENTER_FLUID_COLD;
	}
	if(y == 1)
	{
		Prop_water(Tci);
		Cpc = c;
	}
	if(y == 2)
	{
		Prop_air(Tci);
		Cpc = c;
	}
	if(y == 3)
	{
		Prop_oil(Tci);
		Cpc = c;
	}
	if(y == 4)
	{
		printf("Enter the Specific Heat of Cold Fluid at Temp %lf 'c (in J/Kg.K): ",Tci);
		scanf("%lf",&Cpc);
	}
	fprintf(fp1,"Cph = %lf J/Kg.K . \n",Cph);
	fprintf(fp2,"Cpc = %lf J/Kg.K . \n",Cpc);
	Ch = Mh*Cph;
	Cc = Mc*Cpc;
	fprintf(fp1,"Ch = %lf J/K.sec . \n",Ch);
	fprintf(fp2,"Cc = %lf J/K.sec . \n",Cc);
	fprintf(fp8,"Mh \t %lf \t Kg/sec . \t",Mh);
	fprintf(fp8,"Cph \t %lf \t J/Kg.K . \n",Cph);
	fprintf(fp8,"Mc \t %lf \t Kg/sec . \t",Mc);
	fprintf(fp8,"Cpc \t %lf \t J/Kg.K . \n",Cpc);
	fprintf(fp8,"Ch \t %lf \t J/K.Sec . \t",Ch);
	fprintf(fp8,"Cc \t %lf \t J/K.Sec . \n",Cc);
	fprintf(fp8,"\n");
	fprintf(fp8,"\n");
	fprintf(fp8,"Epsilon \t Tho('c) \t Tco('c) \t Tcommin('c) \t Tcommout('c) \t LMTD \n");
	for(n = 0.1; n < 1; n = n+0.01)
	{
		Tho = CALC_THO(n,Cc,Ch,Thi,Tci);
		Tco = CALC_TCO(n,Cc,Ch,Thi,Tci);
		Tcommin = Thi - Tco;
		Tcommout = Tho - Tci;
		LMTD = CALC_LMTD(Tcommin,Tcommout);
		fprintf(fp8,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",n,Tho,Tco,Tcommin,Tcommout,LMTD);
	}
	fclose(fp8);
	Tho = CALC_THO(Epsilon,Cc,Ch,Thi,Tci);
	Tco = CALC_TCO(Epsilon,Cc,Ch,Thi,Tci);
	fprintf(fp1,"Epsilon = %lf \n",Epsilon);
	fprintf(fp2,"Epsilon = %lf \n",Epsilon);
	fprintf(fp1,"Hot fluid Final Temperature (Tho) = %lf 'c . \n",Tho);
	fprintf(fp2,"Cold fluid Final Temperature (Tco) = %lf 'c . \n",Tco);
	Tbh = (Thi + Tho)/2;
	Tbc = (Tci + Tco)/2;
	fprintf(fp1,"Bulk Mean Temperature for Hot Fluid (Tbh) = %lf 'c . \n",Tbh);
	fprintf(fp2,"Bulk Mean Temperature for Cold Fluid (Tbc) = %lf 'c . \n",Tbc);
	if(y == 1)
	{
		Prop_water(Tbc);
		Cpc = c;
		k_c = k;
		rho_c = r;
		vis_c = v;		
	}
	if(y == 2)
	{
		Prop_air(Tbc);
		Cpc = c;
		k_c = k;
		rho_c = r;
		vis_c = v;
	}
	if(y == 3)
	{
		Prop_oil(Tbc);
		Cpc = c;
		k_c = k;
		rho_c = r;
		vis_c = v;
	}
	if(y == 4)
	{
		printf("Enter the Specific Heat of Cold Fluid at Temp %lf 'c (in J/Kg.K) : ",Tbc);
		scanf("%lf",&Cpc);
		printf("Enter the Thermal Conductivity of Cold Fluid at Temp %lf 'c (in W/m.K) : ",Tbc);
		scanf("%lf",&k_c);
		printf("Enter the Density of Cold Fluid at Temp %lf 'c (in Kg/m3) : ",Tbc);
		scanf("%lf",&rho_c);
		printf("Enter the Viscosity of Cold Fluid at Temp %lf 'c (in Pa.s) : ",Tbc);
		scanf("%lf",&vis_c); 
	}
	fprintf(fp2,"Specific Heat of Cold Fluid (Cpc) = %lf J/Kg.K . \n",Cpc);
	fprintf(fp2,"Thermal Conductivity of Cold Fluid (Kc) = %lf W/m.K . \n",k_c);
	fprintf(fp2,"Density of Cold Fluid (rhoc) = %lf Kg/m3 . \n",rho_c);
	fprintf(fp2,"Viscosity of Cold Fluid (visc) = %lf Pa.s . \n",vis_c);
	if(x == 1)
	{
		Prop_water(Tbh);
		Cph = c;
		k_h = k;
		rho_h = r;
		vis_h = v;
	}
	if(x == 2)
	{
		Prop_air(Tbh);
		Cph = c;
		k_h = k;
		rho_h = r;
		vis_h = v;
	}
	if(x == 3)
	{
		Prop_oil(Tbh);
		Cph = c;
		k_h = k;
		rho_h = r;
		vis_h = v;
	}
	if(x == 4)
	{
		printf("Enter the Specific Heat of Hot Fluid at Temp %lf 'c (in J/Kg.K) : ",Tbh);
		scanf("%lf",&Cph);
		printf("Enter the Thermal Conductivity of Hot Fluid at Temp %lf 'c (in W/m.K) : ",Tbh);
		scanf("%lf",&k_h);
		printf("Enter the Density of Hot Fluid at Temp %lf 'c (in Kg/m3) : ",Tbh);
		scanf("%lf",&rho_h);
		printf("Enter the Viscosity of Hot Fluid at Temp %lf 'c (in Pa.s) : ",Tbh);
		scanf("%lf",&vis_h); 
	}
	fprintf(fp1,"Specific Heat of Hot Fluid (Cph) = %lf J/Kg.K . \n",Cph);
	fprintf(fp1,"Thermal Conductivity of Hot Fluid (Kh) = %lf W/m.K . \n",k_h);
	fprintf(fp1,"Density of Hot Fluid (rhoh) = %lf Kg/m3 . \n",rho_h);
	fprintf(fp1,"Viscosity of Hot Fluid (vish) = %lf Pa.s . \n",vis_h);
	Ch = Mh*Cph;
	Cc = Mc*Cpc;
	fprintf(fp1,"Ch = %lf J/K.sec . \n",Ch);
	fprintf(fp2,"Cc = %lf J/K.sec . \n",Cc);
	ENTER_MATERIAL_SELECTION:
	printf("Choose the Material of Heat Exchanger from below options : \n");
	printf("1. Plain carbon steel \n");
	printf("2. Stainless steel \n");
	printf("3. Aluminum \n");
	printf("4. Copper \n");
	printf("5. 90-10 cupronickel \n");
	printf("6. 70-30 cupronickel \n");
	printf("7. Titanium \n");
	printf("8. Other \n");
	printf("Enter your selected option : ");
	scanf("%d",&z);
	if(z < 1 || z > 8)
	{
		printf("Invalid input of Material. \n");
		goto ENTER_MATERIAL_SELECTION;
	}
	if(z > 0 && z < 8)
	{
		Vmax = MATERIAL_VELOCITY(z);
		k_m = MATERIAL_THERMALCONDUCTIVITY(z);
	}
	if(z == 8)
	{
		printf("Enter the Maximum Allowable Velocity of Material (in m/s) : ");
		scanf("%lf",&Vmax);
		printf("Enter the Thermal Conductivity of Material (in W/mK): ");
		scanf("%lf",&k_m);
	}
	fprintf(fp1,"Maximum Velocity of Material (Vmax) = %lf m/s . \n",Vmax);
	fprintf(fp1,"Thermal Conductivity of Material (Km) = %lf W/m.K . \n",k_m);
	fprintf(fp2,"Maximum Velocity of Material (Vmax) = %lf m/s . \n",Vmax);
	fprintf(fp2,"Thermal Conductivity of Material (Km) = %lf W/m.K . \n",k_m);
	if(Ch < Cc)
	{
		fprintf(fp3,"In Tube side of Heat Exchanger There will be Hot Fluid. \n");
		fprintf(fp3,"Mh = %lf Kg/s . \n",Mh);
		fprintf(fp3,"Specific Heat = %lf J/Kg.K .\n",Cph);
		fprintf(fp3,"Thermal conductivity = %lf W/m.K . \n",k_h);
		fprintf(fp3,"density = %lf Kg/m3 .\n",rho_h);
		fprintf(fp3,"viscosity = %lf Pa.s .\n",vis_h);
		fprintf(fp3,"Effectiveness = %lf \n",Epsilon);
		fprintf(fp3,"Initial Temperature = %lf 'c . \n",Thi);
		fprintf(fp3,"Final Temperature = %lf 'c . \n",Tho);
		fprintf(fp4,"In Annulus side of Heat Exchanger There will be Cold Fluid. \n");
		fprintf(fp4,"Mass flow rate = %lf Kg/s . \n",Mc);
		fprintf(fp4,"Specific Heat = %lf J/Kg.K . \n",Cpc);
		fprintf(fp4,"Thermal conductivity = %lf W/m.K . \n",k_c);
		fprintf(fp4,"density = %lf Kg/m3 . \n",rho_c);
		fprintf(fp4,"viscosity = %lf Pa.s . \n",vis_c);
		fprintf(fp4,"Effectiveness = %lf . \n",Epsilon);
		fprintf(fp4,"Initial Temperature = %lf 'c . \n",Tci);
		fprintf(fp4,"Final Temperature = %lf 'c . \n",Tco);
		fp5 = fopen("Tube_side_selection_Option.txt","w");
		printf("Hot Fluid is flowing through Tube Side. \n");
		fprintf(fp1,"Hot Fluid is flowing through Tube Side. \n");
		fprintf(fp2,"Hot Fluid is flowing through Tube Side. \n");
		ENTER_SELECTION_INNERDIA_HOT:
		printf("For your Selected Heat exchanger Material this are options for flow your Hot Fluid \n");
		fprintf(fp5,"For your Selected Heat exchanger Material this are options for flow your Hot Fluid \n");
		for(i=0;i<75;i++)
		{
			di = TEMA_TUBESIDE_INNERDIA(i);
			Aci = PI*di*di*pow(10,-6)/4;
			Vmi = Mh/(rho_h*Aci);
			if(vmi < Vmax)
			{
				d0 = TEMA_TUBESIDE_OUTERDIA(i);
				ti = (d0 - di)/2;
				printf("%d. Innerdia = %lf mm , Outer Dia = %lf mm & Thickness = %lf mm . \n",i,di,d0,ti);
				fprintf(fp5,"%d. Innerdia = %lf mm, Outer Dia = %lf mm & Thickness = %lf mm .\n",i,di,d0,ti);
			}
		}
		printf("Choose any one size for Hot Fluid Flow : ");
		fprintf(fp5,"Choose any one size for Hot Fluid Flow : ");
		fclose(fp5);
		scanf("%d",&k1);
		if(k1 < 0 || k1 > 75)
		{
			printf("Invalid Selected Option. \n");
			goto ENTER_SELECTION_INNERDIA_HOT;
		}
		di = TEMA_TUBESIDE_INNERDIA(k1);
		d0 = TEMA_TUBESIDE_OUTERDIA(k1);
		Aci = PI*di*di*pow(10,-6)/4;
		Vmi = Mh/(rho_h*Aci);
		if(Vmi > Vmax)
		{
			printf("Invalid Selected Option. \n");
			goto ENTER_SELECTION_INNERDIA_HOT;
		}
		fprintf(fp1,"di = %lf mm . \n",di);
		fprintf(fp1,"do = %lf mm . \n",d0);
		fprintf(fp2,"di = %lf mm . \n",di);
		fprintf(fp2,"do = %lf mm . \n",d0);
		fprintf(fp3,"Inner dia of tube = %lf mm . \n",di);
		fprintf(fp3,"Outer dia of tube = %lf mm . \n",d0);
		fprintf(fp3,"Cross-sectional Area of Tube = %lf m2 . \n",Aci);
		fprintf(fp3,"Velocity of Flow in tube = %lf m/s . \n",Vmi);
		fprintf(fp1,"Tube Area (Ai) = %lf m2 . \n",Aci);
		fprintf(fp1,"Velocity in Tube side (Vi) = %lf m/s . \n",Vmi);
		fprintf(fp2,"Tube Area (Ai) = %lf m2 . \n",Aci);
		fprintf(fp2,"Velocity in Tube side (Vi) = %lf m/s . \n",Vmi);
		D1 = di + d0;
		fprintf(fp4,"Minimum Dh = %lf mm . \n",D1);
		fprintf(fp1,"Minimum Annulus Dia requirement = %lf mm . \n",D1);
		fprintf(fp2,"Minimum Annulus Dia requirement = %lf mm . \n",D1);
		ENTER_ANNULARSIDE_DIA_COLD:
		printf("As per your selected Inside Tube For Annular Pipe Options are as below : \n");
		fprintf(fp6,"As per your selected Inside Tube For Annular Pipe Options are as below : \n");
		for(j = 0;j < 75;j++)
		{
			Di = TEMA_ANNULARSIDE_INNERDIA(j);
			if(Di == D1 || Di > D1)
			{
				D0 = TEMA_ANNULARSIDE_OUTERDIA(j);
				t0 = (D0 - Di)/2;
				printf("%d . Inner Dia. = %lf , Outer Dia = %lf & Thickness = %lf \n",j,Di,D0,t0);
				fprintf(fp6,"%d . Inner Dia. = %lf mm, Outer Dia = %lf mm & Thickness = %lf mm .\n",j,Di,D0,t0);
			}
		}
		printf("Choose any one size for cold Fluid Flow as minimum innerdia required is %lf : ",D1);
		fprintf(fp6,"Choose any one size for cold Fluid Flow as minimum innerdia required is %lf : ",D1);
		scanf("%d",&z1);
		fclose(fp6);
		if(z1 < 1 || z1 > 75)
		{
			printf("Invalid selected option : \n");
			goto ENTER_ANNULARSIDE_DIA_COLD;
		}
		Di = TEMA_ANNULARSIDE_INNERDIA(z1);
		D0 = TEMA_ANNULARSIDE_OUTERDIA(z1);
		if(Di < D1)
		{
			printf("Invalid selected option : \n");
			goto ENTER_ANNULARSIDE_DIA_COLD;
		}
		fprintf(fp1,"Di = %lf mm . \n",Di);
		fprintf(fp1,"Do = %lf mm . \n",D0);
		fprintf(fp2,"Di = %lf mm . \n",Di);
		fprintf(fp2,"Do = %lf mm . \n",D0);
		fprintf(fp4,"From Tema Table Di = %lf mm .\n",Di);
		fprintf(fp4,"Do = %lf mm .\n",D0);
		Dh = Di - d0;
		fprintf(fp4,"Updated Dh = %lf mm .\n",Dh);
		b1 = Di/1000;
		b2 = d0/1000;
		De = ((pow(b1,2) - pow(b2,2))/b2);
		Aco = ((Di*Di*pow(10,-6)) - (d0*d0*pow(10,-6)))*PI/4;
		Vmo = Mc/(rho_c*Aco);
		Rei = CALC_RE(rho_h,Vmi,di,vis_h);
		Reo = CALC_RE(rho_c,Vmo,Dh,vis_c);
		Pri = CALC_PR(vis_h,Cph,k_h);
		Pro = CALC_PR(vis_c,Cpc,k_c);
		Nui = CALC_NU(Rei,Pri);
		Nuo = CALC_NU(Reo,Pro);
		hi = CALC_H(Nui,k_h,di);
		ho = CALC_H(Nuo,k_c,De);
		Rft = 0.000176;
		fprintf(fp4,"Cross-sectional Area = %lf m2 .\n",Aco);
		fprintf(fp4,"Re = %lf \n",Reo);
		fprintf(fp3,"Re = %lf \n",Rei);
		fprintf(fp3,"Pr = %lf \n",Pri);
		fprintf(fp4,"Pr = %lf \n",Pro);
		fprintf(fp4,"Nu = %lf \n",Nuo);
		fprintf(fp3,"Nu = %lf \n",Nui);
		fprintf(fp3,"hi = %lf W/(K. m2) . \n",hi);
		fprintf(fp4,"De = %lf mm . \n",De);
		fprintf(fp4,"ho = %lf W/(K. m2) . \n",ho);
		Uc = CALC_UC(d0,di,hi,k_m,ho);
		Uf = CALC_UF(Uc,Rft);
		CF = CALC_CF(Uf,Uc);
		OD = CALC_OD(CF);
		fprintf(fp7,"Uc = %lf W/(K . m2) .\n",Uc);
		fprintf(fp7,"Assume Rft = %lf (K . m2)/W . \n",Rft);
		fprintf(fp7,"Uf = %lf W/(K . m2) . \n",Uf);
		fprintf(fp7,"CF = %lf \n",CF);
		fprintf(fp7,"OD = %lf \n",OD);
		Q = Ch*(Thi - Tho);
		Tcommin = Thi - Tco;
		Tcommout = Tho - Tci;
		LMTD = CALC_LMTD(Tcommin,Tcommout);
		fprintf(fp7,"Q = %lf W . \n",Q);
		fprintf(fp7,"Tcommin = %lf 'c . \n",Tcommin);
		fprintf(fp7,"Tcommout = %lf 'c . \n",Tcommout);
		fprintf(fp7,"LMTD = %lf \n",LMTD);
		fprintf(fp1,"Dh = %lf mm . \n",Dh);
		fprintf(fp1,"De = %lf mm . \n",De);
		fprintf(fp1,"Annulus Tube Area (Aco) = %lf m2 . \n",Aco);
		fprintf(fp1,"Velocity of Fluid in Annulus Tube (Vo) = %lf m/s. \n",Vmo);
		fprintf(fp1,"Tube side Re = %lf \n",Rei);
		fprintf(fp1,"Annulus side Re = %lf \n",Reo);
		fprintf(fp1,"Tube side Pr = %lf \n",Pri);
		fprintf(fp1,"Annulus side Pr = %lf \n",Pro);
		fprintf(fp1,"Tube side Nu = %lf \n",Nui);
		fprintf(fp1,"Annulus side Nu = %lf \n",Nuo);
		fprintf(fp1,"hi = %lf W/(K . m2) . \n",hi);
		fprintf(fp1,"ho = %lf W/(K . m2) . \n",ho);
		fprintf(fp1,"Uc = %lf W/(K . m2) . \n",Uc);
		fprintf(fp1,"Uf = %lf W/(K . m2) . \n",Uf);
		fprintf(fp1,"CF = %lf \n",CF);
		fprintf(fp1,"OD = %lf % . \n",OD);
		fprintf(fp1,"Q = %lf W . \n",Q);
		fprintf(fp1,"Tcommin = %lf 'c . \n",Tcommin);
		fprintf(fp1,"Tcommout = %lf 'c . \n",Tcommout);
		fprintf(fp2,"LMTD = %lf \n",LMTD);
		fprintf(fp2,"Dh = %lf mm . \n",Dh);
		fprintf(fp2,"De = %lf mm . \n",De);
		fprintf(fp2,"Annulus Tube Area (Aco) = %lf m2 . \n",Aco);
		fprintf(fp2,"Velocity of Fluid in Annulus Tube (Vo) = %lf m/s . \n",Vmo);
		fprintf(fp2,"Tube side Re = %lf \n",Rei);
		fprintf(fp2,"Annulus side Re = %lf \n",Reo);
		fprintf(fp2,"Tube side Pr = %lf \n",Pri);
		fprintf(fp2,"Annulus side Pr = %lf \n",Pro);
		fprintf(fp2,"Tube side Nu = %lf \n",Nui);
		fprintf(fp2,"Annulus side Nu = %lf \n",Nuo);
		fprintf(fp2,"hi = %lf \n",hi);
		fprintf(fp2,"ho = %lf \n",ho);
		fprintf(fp2,"Uc = %lf \n",Uc);
		fprintf(fp2,"Uf = %lf \n",Uf);
		fprintf(fp2,"CF = %lf \n",CF);
		fprintf(fp2,"OD = %lf % \n",OD);
		fprintf(fp2,"Q = %lf W . \n",Q);
		fprintf(fp2,"Tcommin = %lf 'c . \n",Tcommin);
		fprintf(fp2,"Tcommout = %lf 'c . '\n",Tcommout);
		fprintf(fp2,"LMTD = %lf \n",LMTD);
		Ao = Q/(Uf*LMTD);
		N = Ao/(2*PI*d0*L);
		N1 = N / 1;
		N2 = N - N1;
		fprintf(fp7,"Ao = %lf m2 . \n",Ao);
		fprintf(fp7,"N = %lf \n",N);
		fprintf(fp1,"Outside Area Ao = %lf m2 . \n",Ao);
		fprintf(fp1,"Number of Tubes N = %lf \n",N);
		fprintf(fp2,"Outside Area Ao = %lf m2 . \n",Ao);
		fprintf(fp2,"Number of Tubes N = %lf \n",N);
		if(N2 > 0.5 || N2 == 0.5)
		{
			N = N1 + 1;
		}
		if(N2 < 0.5 || N2 == 0.0)
		{
			N = N1;
		}
		fprintf(fp7,"N ~= %lf \n",N);
		fprintf(fp1,"N ~= %lf \n",N);
		fprintf(fp2,"N ~= %lf \n",N);
		printf("Enter the dimensional roughness of Material of Tubein (in mm) : ");
		scanf("%lf",&ki);
		printf("Enter the dimensional roughness of Material of Annulus (in mm) : ");
		scanf("%lf",&ko);
		fi = calc_f(Rei,ki,di);
		f0 = calc_f(Reo,ko,Di);
		Pin = calc_pressuredrop_Tubein(fi,N,L,di,Vmi,rho_h);
		Pout = calc_pressuredrop_Annulus(f0,N,L,Dh,Vmo,rho_c,h);
		P = Pin + Pout;
		fprintf(fp7,"fi = %lf \n",fi);
		fprintf(fp7,"fo = %lf \n",f0);
		fprintf(fp7,"Pressure drop in Tubeside Pin = %lf Pa . \n",Pin);
		fprintf(fp7,"Pressure drop in Annulus Tube Pout = %lf Pa . \n",Pout);
		fprintf(fp7,"Total Pressure drop P = %lf Pa . \n",P);
		fprintf(fp1,"fi = %lf \n",fi);
		fprintf(fp1,"fo = %lf \n",f0);
		fprintf(fp1,"Pressure drop in Tubeside Pin = %lf Pa . \n",Pin);
		fprintf(fp1,"Pressure drop in Annulus Tube Pout = %lf Pa . \n",Pout);
		fprintf(fp1,"Total Pressure drop P = %lf Pa . \n",P);
		fprintf(fp2,"fi = %lf \n",fi);
		fprintf(fp2,"fo = %lf \n",f0);
		fprintf(fp2,"Pressure drop in Tubeside Pin = %lf Pa . \n",Pin);
		fprintf(fp2,"Pressure drop in Annulus Tube Pout = %lf Pa . \n",Pout);
		fprintf(fp2,"Total Pressure drop P = %lf Pa . \n",P);
		fprintf(fp9,"Q(W) \t Area(m2) \t Length(m) \t Number_of_Tube \t fi \t fo \t Pin(Pa) \t Pout(Pa) \t Ptotal(Pa) \n");
		double m;
		for(m = 0.5; m < 10; m=m+0.5)
		{
			Q = Ch*(Thi - Tho);
			Ao = Q/(Uf*LMTD);
			N = Ao/(2*PI*d0*m);
			N1 = N / 1;
			N2 = N - N1;	
			if(N2 > 0.5 || N2 == 0.5)
			{
				N = N1 + 1;
			}
			if(N2 < 0.5 || N2 == 0.0)
			{
				N = N1;
			}
			fi = calc_f(Rei,ki,di);
			f0 = calc_f(Reo,ko,Di);
			printf("Enter the R/D ratio for Annulus Side pressure drop calculation for 90' bend : ");
			scanf("%lf",&h);
			Pin = calc_pressuredrop_Tubein(fi,N,m,di,Vmi,rho_h);
			Pout = calc_pressuredrop_Annulus(f0,N,m,Dh,Vmo,rho_c,h);
			P = Pin + Pout;
			fprintf(fp9,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",Q,Ao,m,N,fi,f0,Pin,Pout,P);
		}
	}
	if(Cc < Ch)
	{
		fprintf(fp3,"In Tube side of Heat Exchanger There will be Cold Fluid. \n");
		fprintf(fp3,"Mass flow rate = %lf Kg/s .\n",Mc);
		fprintf(fp3,"Specific Heat = %lf J/Kg.K .\n",Cpc);
		fprintf(fp3,"Thermal conductivity = %lf W/m.K .\n",k_c);
		fprintf(fp3,"density = %lf Kg/m3 .\n",rho_c);
		fprintf(fp3,"viscosity = %lf Pa.s .\n",vis_c);
		fprintf(fp3,"Effectiveness = %lf \n",Epsilon);
		fprintf(fp3,"Initial Temperature = %lf 'c \n",Tci);
		fprintf(fp3,"Final Temperature = %lf 'c \n",Tco);
		fprintf(fp4,"In Annulus side of Heat Exchanger There will be Hot Fluid. \n");
		fprintf(fp4,"Mass flow rate = %lf Kg/s \n",Mh);
		fprintf(fp4,"Specific Heat = %lf J/Kg.K .\n",Cph);
		fprintf(fp4,"Thermal conductivity = %lf W/m.K .\n",k_h);
		fprintf(fp4,"density = %lf Kg/m3 .\n",rho_h);
		fprintf(fp4,"viscosity = %lf Pa.s .\n",vis_h);
		fprintf(fp4,"Effectiveness = %lf \n",Epsilon);
		fprintf(fp4,"Initial Temperature = %lf 'c \n",Thi);
		fprintf(fp4,"Final Temperature = %lf 'c \n",Tho);
		ENTER_SELECTION_INNERDIA_COLD:
		printf("For your Selected Heat exchanger Material this are options for flow your Cols Fluid \n");
		fprintf(fp5,"For your Selected Heat exchanger Material this are options for flow your Cols Fluid \n");
		for(i=0;i<75;i++)
		{
			di = TEMA_TUBESIDE_INNERDIA(i);
			Aci = PI*di*di*pow(10,-6)/4;
			Vmi = Mc/(rho_c*Aci);
			if(vmi < Vmax)
			{
				d0 = TEMA_TUBESIDE_OUTERDIA(i);
				ti = (d0 - di)/2;
				printf("%d . Innerdia = %lf , Outer Dia = %lf & Thickness = %lf \n",i,di,d0,ti);
				fprintf(fp5,"%d . Innerdia = %lf , Outer Dia = %lf & Thickness = %lf \n",i,di,d0,ti);
			}
		}
		printf("Choose any one size for Hot Fluid Flow : ");
		fprintf(fp5,"Choose any one size for Hot Fluid Flow : ");
		fclose(fp5);
		scanf("%d",&k1);
		if(k1 < 1 || k1 > 75)
		{
			printf("Invalid Selected Option. \n");
			goto ENTER_SELECTION_INNERDIA_COLD;
		}
		di = TEMA_TUBESIDE_INNERDIA(k1);
		d0 = TEMA_TUBESIDE_OUTERDIA(k1);
		Aci = PI*di*di*pow(10,-6)/4;
		Vmi = Mc/(rho_c*Aci);
		if(Vmi > Vmax)
		{
			printf("Invalid Selected Option. \n");
			goto ENTER_SELECTION_INNERDIA_COLD;
		}
		fprintf(fp3,"Inner dia of tube = %lf mm . \n",di);
		fprintf(fp3,"Outer dia of tube = %lf mm . \n",d0);
		fprintf(fp3,"Cross-sectional Area of Tube = %lf m2 . \n",Aci);
		fprintf(fp3,"Velocity of Flow in tube = %lf m/s . \n",Vmi);
		fprintf(fp1,"di = %lf mm . \n",di);
		fprintf(fp1,"do = %lf mm . \n",d0);
		fprintf(fp2,"di = %lf mm . \n",di);
		fprintf(fp2,"do = %lf mm . \n",d0);
		fprintf(fp1,"Tube Area (Ai) = %lf m2 . \n",Aci);
		fprintf(fp1,"Velocity in Tube side (Vi) = %lf m/s . \n",Vmi);
		fprintf(fp2,"Tube Area (Ai) = %lf m2 . \n",Aci);
		fprintf(fp2,"Velocity in Tube side (Vi) = %lf m/s . \n",Vmi);
		D1 = di + d0;
		fprintf(fp4,"Minimum Dh = %lf mm . \n",D1);
		fprintf(fp1,"Minimum Annulus Dia requirement = %lf mm . \n",D1);
		fprintf(fp2,"Minimum Annulus Dia requirement = %lf mm . \n",D1);
		ENTER_ANNULARSIDE_DIA_HOT:
		printf("As per your selected Inside Tube For Annular Pipe Options are as below : \n");
		fprintf(fp6,"As per your selected Inside Tube For Annular Pipe Options are as below : \n");
		for(j = 0;j < 75;j++)
		{
			Di = TEMA_ANNULARSIDE_INNERDIA(j);
			if(Di == D1 || Di > D1)
			{
				D0 = TEMA_ANNULARSIDE_OUTERDIA(j);
				t0 = (D0 - Di)/2;
				printf("%d . Inner Dia. = %lf mm, Outer Dia = %lf mm & Thickness = %lf mm . \n",j,Di,D0,t0);
				fprintf(fp6,"%d . Inner Dia. = %lf mm, Outer Dia = %lf mm & Thickness = %lf mm . \n",j,Di,D0,t0);
			}
		}
		printf("Choose any one size for cold Fluid Flow as minimum innerdia required is %lf : ",D1);
		fprintf(fp6,"Choose any one size for cold Fluid Flow as minimum innerdia required is %lf : ",D1);
		scanf("%d",&z1);
		if(z1 < 1 || z1 > 75)
		{
			printf("Invalid selected option : \n");
			goto ENTER_ANNULARSIDE_DIA_HOT;
		}
		Di = TEMA_ANNULARSIDE_INNERDIA(z1);
		D0 = TEMA_ANNULARSIDE_OUTERDIA(z1);
		if(Di < D1)
		{
			printf("Invalid selected option : \n");
			goto ENTER_ANNULARSIDE_DIA_HOT;
		}
		fprintf(fp4,"From Tema Table Di = %lf mm . \n",Di);
		fprintf(fp4,"Do = %lf mm . \n",D0);
		Dh = Di - d0;
		fprintf(fp4,"Updated Dh = %lf mm . \n",Dh);
		fprintf(fp1,"Di = %lf mm . \n",Di);
		fprintf(fp1,"Do = %lf mm . \n",D0);
		fprintf(fp2,"Di = %lf mm . \n",Di);
		fprintf(fp2,"Do = %lf mm . \n",D0);
		b1 = Di/1000;
		b2 = d0/1000;
		De = ((pow(b1,2) - pow(b2,2))/b2);
		De = De*1000;
		Aco = ((Di*Di*pow(10,-6)) - (d0*d0*pow(10,-6)))*PI/4;
		Vmo = Mh/(rho_h*Aco);
		Rei = CALC_RE(rho_c,Vmi,di,vis_c);
		Reo = CALC_RE(rho_h,Vmo,Dh,vis_h);
		Pri = CALC_PR(vis_c,Cpc,k_c);
		Pro = CALC_PR(vis_h,Cph,k_h);
		Nui = CALC_NU(Rei,Pri);
		Nuo = CALC_NU(Reo,Pro);
		hi = CALC_H(Nui,k_c,di);
		ho = CALC_H(Nuo,k_h,De);
		fprintf(fp4,"Cross-sectional Area = %lf m2 . \n",Aco);
		fprintf(fp4,"Re = %lf \n",Reo);
		fprintf(fp3,"Re = %lf \n",Rei);
		fprintf(fp3,"Pr = %lf \n",Pri);
		fprintf(fp4,"Pr = %lf \n",Pro);
		fprintf(fp4,"Nu = %lf \n",Nuo);
		fprintf(fp3,"Nu = %lf \n",Nui);
		fprintf(fp3,"hi = %lf W/(K . m2) . \n",hi);
		fprintf(fp4,"De = %lf mm .\n",De);
		fprintf(fp4,"ho = %lf W/(K . m2) . \n",ho);
		Rft = 0.000176;
		Uc = CALC_UC(d0,di,hi,k_m,ho);
		Uf = CALC_UF(Uc,Rft);
		CF = CALC_CF(Uf,Uc);
		OD = CALC_OD(CF);
		Q = Ch*(Thi - Tho);
		Tcommin = Thi - Tco;
		Tcommout = Tho - Tci;
		LMTD = CALC_LMTD(Tcommin,Tcommout);
		fprintf(fp7,"Q = %lf W . \n",Q);
		fprintf(fp7,"Tcommin = %lf 'c . \n",Tcommin);
		fprintf(fp7,"Tcommout = %lf 'c . \n",Tcommout);
		fprintf(fp7,"LMTD = %lf \n",LMTD);
		fprintf(fp1,"Dh = %lf mm . \n",Dh);
		fprintf(fp1,"De = %lf mm . \n",De);
		fprintf(fp1,"Annulus Tube Area (Aco) = %lf m2 . \n",Aco);
		fprintf(fp1,"Velocity of Fluid in Annulus Tube (Vo) = %lf m/s . \n",Vmo);
		fprintf(fp1,"Tube side Re = %lf \n",Rei);
		fprintf(fp1,"Annulus side Re = %lf \n",Reo);
		fprintf(fp1,"Tube side Pr = %lf \n",Pri);
		fprintf(fp1,"Annulus side Pr = %lf \n",Pro);
		fprintf(fp1,"Tube side Nu = %lf \n",Nui);
		fprintf(fp1,"Annulus side Nu = %lf \n",Nuo);
		fprintf(fp1,"hi = %lf W/(K . m2) . \n",hi);
		fprintf(fp1,"ho = %lf W/(K . m2) . \n",ho);
		fprintf(fp1,"Uc = %lf W/(K . m2) . \n",Uc);
		fprintf(fp1,"Uf = %lf W/(K . m2) . \n",Uf);
		fprintf(fp1,"CF = %lf . \n",CF);
		fprintf(fp1,"OD = %lf % . \n",OD);
		fprintf(fp1,"Q = %lf W . \n",Q);
		fprintf(fp1,"Tcommin = %lf 'c . \n",Tcommin);
		fprintf(fp1,"Tcommout = %lf 'c . \n",Tcommout);
		fprintf(fp2,"LMTD = %lf \n",LMTD);
		fprintf(fp2,"Dh = %lf mm . \n",Dh);
		fprintf(fp2,"De = %lf mm . \n",De);
		fprintf(fp2,"Annulus Tube Area (Aco) = %lf m2 . \n",Aco);
		fprintf(fp2,"Velocity of Fluid in Annulus Tube (Vo) = %lf m/s . \n",Vmo);
		fprintf(fp2,"Tube side Re = %lf \n",Rei);
		fprintf(fp2,"Annulus side Re = %lf \n",Reo);
		fprintf(fp2,"Tube side Pr = %lf \n",Pri);
		fprintf(fp2,"Annulus side Pr = %lf \n",Pro);
		fprintf(fp2,"Tube side Nu = %lf \n",Nui);
		fprintf(fp2,"Annulus side Nu = %lf \n",Nuo);
		fprintf(fp2,"hi = %lf W/(K . m2) . \n",hi);
		fprintf(fp2,"ho = %lf W/(K . m2) . \n",ho);
		fprintf(fp2,"Uc = %lf W/(K . m2) . \n",Uc);
		fprintf(fp2,"Uf = %lf W/(K . m2) . \n",Uf);
		fprintf(fp2,"CF = %lf \n",CF);
		fprintf(fp2,"OD = %lf % . \n",OD);
		fprintf(fp2,"Q = %lf W . \n",Q);
		fprintf(fp2,"Tcommin = %lf 'c . \n",Tcommin);
		fprintf(fp2,"Tcommout = %lf 'c . \n",Tcommout);
		fprintf(fp2,"LMTD = %lf \n",LMTD);
		Ao = Q/(Uf*LMTD);
		N = Ao/(2*PI*d0*L);
		N1 = N / 1;
		N2 = N - N1;
		fprintf(fp7,"Ao = %lf m2 . \n",Ao);
		fprintf(fp7,"N = %lf . \n",N);
		fprintf(fp1,"Outside Area Ao = %lf m2 . \n",Ao);
		fprintf(fp1,"Number of Tubes N = %lf . \n",N);
		fprintf(fp2,"Outside Area Ao = %lf m2 . \n",Ao);
		fprintf(fp2,"Number of Tubes N = %lf \n",N);
		if(N2 > 0.5 || N2 == 0.5)
		{
			N = N1 + 1;
		}
		if(N2 < 0.5 || N2 == 0.0)
		{
			N = N1;
		}
		fprintf(fp7,"N ~= %lf \n",N);
		fprintf(fp1,"N ~= %lf \n",N);
		fprintf(fp2,"N ~= %lf \n",N);
		printf("Enter the dimensional roughness of Material of Tubein (in mm) : ");
		scanf("%lf",&ki);
		printf("Enter the dimensional roughness of Material of Annulus (in mm) : ");
		scanf("%lf",&ko);
		fi = calc_f(Rei,ki,di);
		f0 = calc_f(Reo,ko,De);
		printf("Enter the R/D ratio for Annulus Side pressure drop calculation for 90' bend : ");
		scanf("%lf",&h);
		Pin = calc_pressuredrop_Tubein(fi,N,L,di,Vmi,rho_c);
		Pout = calc_pressuredrop_Annulus(f0,N,L,Dh,Vmo,rho_h,h);
		P = Pin + Pout;
		fprintf(fp7,"fi = %lf . \n",fi);
		fprintf(fp7,"fo = %lf .\n",f0);
		fprintf(fp7,"Pressure drop in Tubeside Pin = %lf Pa .\n",Pin);
		fprintf(fp7,"Pressure drop in Annulus Tube Pout = %lf Pa . \n",Pout);
		fprintf(fp7,"Total Pressure drop P = %lf Pa . \n",P);
		fprintf(fp1,"fi = %lf \n",fi);
		fprintf(fp1,"fo = %lf \n",f0);
		fprintf(fp1,"Pressure drop in Tubeside Pin = %lf Pa . \n",Pin);
		fprintf(fp1,"Pressure drop in Annulus Tube Pout = %lf Pa . \n",Pout);
		fprintf(fp1,"Total Pressure drop P = %lf Pa . \n",P);
		fprintf(fp2,"fi = %lf \n",fi);
		fprintf(fp2,"fo = %lf \n",f0);
		fprintf(fp2,"Pressure drop in Tubeside Pin = %lf Pa . \n",Pin);
		fprintf(fp2,"Pressure drop in Annulus Tube Pout = %lf Pa . \n",Pout);
		fprintf(fp2,"Total Pressure drop P = %lf Pa . \n",P);
		fprintf(fp9,"Q(W) \t Area(m2) \t Length(m) \t Number_of_Tube \t fi \t fo \t Pin(Pa) \t Pout(Pa) \t Ptotal(Pa) \n");
		double m;
		for(m = 0.5; m < 10; m=m+0.5)
		{
			Q = Ch*(Thi - Tho);
			Ao = Q/(Uf*LMTD);
			N = Ao/(2*PI*d0*m);
			N1 = N / 1;
			N2 = N - N1;	
			if(N2 > 0.5 || N2 == 0.5)
			{
				N = N1 + 1;
			}
			if(N2 < 0.5 || N2 == 0.0)
			{
				N = N1;
			}
			fi = calc_f(Rei,ki,di);
			f0 = calc_f(Reo,ko,Di);
			Pin = calc_pressuredrop_Tubein(fi,N,m,di,Vmi,rho_h);
			Pout = calc_pressuredrop_Annulus(f0,N,m,Dh,Vmo,rho_c,h);
			P = Pin + Pout;
			fprintf(fp9,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",Q,Ao,m,N,fi,f0,Pin,Pout,P);
		}
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp7);
	fclose(fp9);
}
double CALC_TCO(double e, double cc, double ch, double th, double tc)
{
	double q,tco,cm;
	if(ch < cc)
	{
		cm = ch;
	}
	if(cc < ch)
	{
		cm = cc;
	}
	q = e*cm*(th - tc);
	tco = tc + (q/cc);
	return tco;
}
double CALC_THO(double e, double cc, double ch, double th, double tc)
{
	double q,tho,cm;
	if(ch < cc)
	{
		cm = ch;
	}
	if(cc < ch)
	{
		cm = cc;
	}
	q = e*cm*(th - tc);
	tho = th - (q/ch);
	return tho;
}
double MATERIAL_VELOCITY(int a)
{
	a5 = fopen(".\\Tube Material\\material_max_velocity.txt","r");
	double vt,b[10],v[10];
	int i;
	for(i = 0;i < 7; i++)
	{
		fscanf(a5,"%lf %lf",&b[i],&v[i]);
		if(a == b[i])
		{
			vt = v[i];
		}
	}
	fclose(a5);
	return vt;
}
double MATERIAL_THERMALCONDUCTIVITY(int c)
{
	a4 = fopen(".\\Tube Material\\Material_Thermal_Conductivity.txt","r");
	double kt,b[10],v[10];
	int i;
	for(i = 0;i < 7; i++)
	{
		fscanf(a4,"%lf %lf",&b[i],&v[i]);
		if(c == b[i])
		{
			kt = v[i];
		}
	}
	fclose(a4);
	return kt;
}
double TEMA_TUBESIDE_INNERDIA(int b)
{
	double di;
	if(b == 1)
	{
		di = 4.93;
	}
	if(b == 2)
	{
		di = 5.23;
	}
	if(b == 3)
	{
		di = 5.44;
	}
	if(b == 4)
	{
		di = 5.54;
	}
	if(b == 5)
	{
		di = 7.04;
	}
	if(b == 6)
	{
		di = 7.75;
	}
	if(b == 7)
	{
		di = 8.11;
	}
	if(b == 8)
	{
		di = 8.41;
	}
	if(b == 9)
	{
		di = 9.40;
	}
	if(b == 10)
	{
		di = 10.21;
	}
	if(b == 11)
	{
		di = 10.92;
	}
	if(b == 12)
	{
		di = 11.28;
	}
	if(b == 13)
	{
		di = 10.34;
	}
	if(b == 14)
	{
		di = 11.05;
	}
	if(b == 15)
	{
		di = 11.66;
	}
	if(b == 16)
	{
		di = 12.22;
	}
	if(b == 17)
	{
		di = 12.58;
	}
	if(b == 18)
	{
		di = 12.93;
	}
	if(b == 19)
	{
		di = 13.39;
	}
	if(b == 20)
	{
		di = 13.75;
	}
	if(b == 21)
	{
		di = 14.10;
	}
	if(b == 22)
	{
		di = 12.24;
	}
	if(b == 23)
	{
		di = 12.95;
	}
	if(b == 24)
	{
		di = 13.51;
	}
	if(b == 25)
	{
		di = 14.22;
	}
	if(b == 26)
	{
		di = 14.83;
	}
	if(b == 27)
	{
		di = 15.39;
	}
	if(b == 28)
	{
		di = 15.75;
	}
	if(b == 29)
	{
		di = 16.10;
	}
	if(b == 30)
	{
		di = 16.56;
	}
	if(b == 31)
	{
		di = 17.27;
	}
	if(b == 32)
	{
		di = 15.42;
	}
	if(b == 33)
	{
		di = 16.13;
	}
	if(b == 34)
	{
		di = 16.69;
	}
	if(b == 35)
	{
		di = 17.40;
	}
	if(b == 36)
	{
		di = 18.01;
	}
	if(b == 37)
	{
		di = 18.57;
	}
	if(b == 38)
	{
		di = 18.93;
	}
	if(b == 39)
	{
		di = 19.28;
	}
	if(b == 40)
	{
		di = 19.74;
	}
	if(b == 41)
	{
		di = 20.45;
	}
	if(b == 42)
	{
		di = 17.02;
	}
	if(b == 43)
	{
		di = 18.59;
	}
	if(b == 44)
	{
		di = 19.30;
	}
	if(b == 45)
	{
		di = 19.86;
	}
	if(b == 46)
	{
		di = 20.57;
	}
	if(b == 47)
	{
		di = 21.18;
	}
	if(b == 48)
	{
		di = 21.74;
	}
	if(b == 49)
	{
		di = 22.10;
	}
	if(b == 50)
	{
		di = 22.91;
	}
	if(b == 51)
	{
		di = 23.62;
	}
	if(b == 52)
	{
		di = 22.61;
	}
	if(b == 53)
	{
		di = 23.37;
	}
	if(b == 54)
	{
		di = 24.94;
	}
	if(b == 55)
	{
		di = 25.65;
	}
	if(b == 56)
	{
		di = 26.21;
	}
	if(b == 57)
	{
		di = 26.92;
	}
	if(b == 58)
	{
		di = 27.53;
	}
	if(b == 59)
	{
		di = 28.45;
	}
	if(b == 60)
	{
		di = 29.26;
	}
	if(b == 61)
	{
		di = 29.97;
	}
	if(b == 62)
	{
		di = 31.29;
	}
	if(b == 63)
	{
		di = 32.56;
	}
	if(b == 64)
	{
		di = 33.88;
	}
	if(b == 65)
	{
		di = 34.80;
	}
	if(b == 66)
	{
		di = 44.70;
	}
	if(b == 67)
	{
		di = 45.26;
	}
	if(b == 68)
	{
		di = 45.97;
	}
	if(b == 69)
	{
		di = 46.58;
	}
	if(b == 70)
	{
		di = 56.69;
	}
	if(b == 71)
	{
		di = 57.96;
	}
	if(b == 72)
	{
		di = 59.28;
	}
	if(b == 73)
	{
		di = 69.39;
	}
	if(b == 74)
	{
		di = 70.66;
	}
	if(b == 75)
	{
		di = 71.98;
	}
	return di;
}
double TEMA_TUBESIDE_OUTERDIA(int i)
{
	double d0;
	if(i == 1 || i == 2 || i == 3 || i == 4)
	{
		d0 = 6.35;
	}
	if(i == 5 || i == 6 || i == 7 || i == 8)
	{
		d0 = 9.53;
	}
	if(i == 9 || i == 10 || i == 11 || i == 12)
	{
		d0 = 12.7;
	}
	if(i == 13 || i == 14 || i == 15 || i == 16 || i == 17 || i == 18 || i == 19 || i == 20 || i == 21)
	{
		d0 = 15.88;
	}
	if(i == 22 || i == 23 || i == 24 || i == 25 || i == 26 || i == 27 || i == 28 || i == 29 || i == 30 || i == 31)
	{
		d0 = 19.05;
	}
	if(i == 32 || i == 33 || i == 34 || i == 35 || i == 36 || i == 37 || i == 38 || i == 39 || i == 40 || i == 41)
	{
		d0 = 22.23;
	}
	if(i == 42 || i == 43 || i == 44 || i == 45 || i == 46 || i == 47 || i == 48 || i == 49 || i == 50 || i == 51)
	{
		d0 = 25.4;
	}
	if(i == 52 || i == 53 || i == 54 || i == 55 || i == 56 || i == 57 || i == 58 || i == 59 || i == 60 || i == 61)
	{
		d0 = 31.75;
	}
	if(i == 62 || i == 63 || i == 64 || i == 65)
	{
		d0 = 38.1;
	}
	if(i == 66 || i == 67 || i == 68 || i == 69)
	{
		d0 = 50.8;
	}
	if(i == 70 || i == 71 || i == 72)
	{
		d0 = 63.5;
	}
	if(i == 73 || i == 74 || i == 75)
	{
		d0 = 76.2;
	}
	return d0;
}
double TEMA_ANNULARSIDE_INNERDIA(int a)
{
	double Di;
	if(a == 1)
	{
		Di = 4.93;
	}
	if(a == 2)
	{
		Di = 5.23;
	}
	if(a == 3)
	{
		Di = 5.44;
	}
	if(a == 4)
	{
		Di = 5.54;
	}
	if(a == 5)
	{
		Di = 7.04;
	}
	if(a == 6)
	{
		Di = 7.75;
	}
	if(a == 7)
	{
		Di = 8.11;
	}
	if(a == 8)
	{
		Di = 8.41;
	}
	if(a == 9)
	{
		Di = 9.40;
	}
	if(a == 10)
	{
		Di = 10.21;
	}
	if(a == 11)
	{
		Di = 10.92;
	}
	if(a == 12)
	{
		Di = 11.28;
	}
	if(a == 13)
	{
		Di = 10.34;
	}
	if(a == 14)
	{
		Di = 11.05;
	}
	if(a == 15)
	{
		Di = 11.66;
	}
	if(a == 16)
	{
		Di = 12.22;
	}
	if(a == 17)
	{
		Di = 12.58;
	}
	if(a == 18)
	{
		Di = 12.93;
	}
	if(a == 19)
	{
		Di = 13.39;
	}
	if(a == 20)
	{
		Di = 13.75;
	}
	if(a == 21)
	{
		Di = 14.10;
	}
	if(a == 22)
	{
		Di = 12.24;
	}
	if(a == 23)
	{
		Di = 12.95;
	}
	if(a == 24)
	{
		Di = 13.51;
	}
	if(a == 25)
	{
		Di = 14.22;
	}
	if(a == 26)
	{
		Di = 14.83;
	}
	if(a == 27)
	{
		Di = 15.39;
	}
	if(a == 28)
	{
		Di = 15.75;
	}
	if(a == 29)
	{
		Di = 16.10;
	}
	if(a == 30)
	{
		Di = 16.56;
	}
	if(a == 31)
	{
		Di = 17.27;
	}
	if(a == 32)
	{
		Di = 15.42;
	}
	if(a == 33)
	{
		Di = 16.13;
	}
	if(a == 34)
	{
		Di = 16.69;
	}
	if(a == 35)
	{
		Di = 17.40;
	}
	if(a == 36)
	{
		Di = 18.01;
	}
	if(a == 37)
	{
		Di = 18.57;
	}
	if(a == 38)
	{
		Di = 18.93;
	}
	if(a == 39)
	{
		Di = 19.28;
	}
	if(a == 40)
	{
		Di = 19.74;
	}
	if(a == 41)
	{
		Di = 20.45;
	}
	if(a == 42)
	{
		Di = 17.02;
	}
	if(a == 43)
	{
		Di = 18.59;
	}
	if(a == 44)
	{
		Di = 19.30;
	}
	if(a == 45)
	{
		Di = 19.86;
	}
	if(a == 46)
	{
		Di = 20.57;
	}
	if(a == 47)
	{
		Di = 21.18;
	}
	if(a == 48)
	{
		Di = 21.74;
	}
	if(a == 49)
	{
		Di = 22.10;
	}
	if(a == 50)
	{
		Di = 22.91;
	}
	if(a == 51)
	{
		Di = 23.62;
	}
	if(a == 52)
	{
		Di = 22.61;
	}
	if(a == 53)
	{
		Di = 23.37;
	}
	if(a == 54)
	{
		Di = 24.94;
	}
	if(a == 55)
	{
		Di = 25.65;
	}
	if(a == 56)
	{
		Di = 26.21;
	}
	if(a == 57)
	{
		Di = 26.92;
	}
	if(a == 58)
	{
		Di = 27.53;
	}
	if(a == 59)
	{
		Di = 28.45;
	}
	if(a == 60)
	{
		Di = 29.26;
	}
	if(a == 61)
	{
		Di = 29.97;
	}
	if(a == 62)
	{
		Di = 31.29;
	}
	if(a == 63)
	{
		Di = 32.56;
	}
	if(a == 64)
	{
		Di = 33.88;
	}
	if(a == 65)
	{
		Di = 34.80;
	}
	if(a == 66)
	{
		Di = 44.70;
	}
	if(a == 67)
	{
		Di = 45.26;
	}
	if(a == 68)
	{
		Di = 45.97;
	}
	if(a == 69)
	{
		Di = 46.58;
	}
	if(a == 70)
	{
		Di = 56.69;
	}
	if(a == 71)
	{
		Di = 57.96;
	}
	if(a == 72)
	{
		Di = 59.28;
	}
	if(a == 73)
	{
		Di = 69.39;
	}
	if(a == 74)
	{
		Di = 70.66;
	}
	if(a == 75)
	{
		Di = 71.98;
	}
	return Di;
}
double TEMA_ANNULARSIDE_OUTERDIA(int c)
{
	double D0;
	if(c == 1 || c == 2 || c == 3 || c == 4)
	{
		D0 = 6.35;
	}
	if(c == 5 || c == 6 || c == 7 || c == 8)
	{
		D0 = 9.53;
	}
	if(c == 9 || c == 10 || c == 11 || c == 12)
	{
		D0 = 12.7;
	}
	if(c == 13 || c == 14 || c == 15 || c == 16 || c == 17 || c == 18 || c == 19 || c == 20 || c == 21)
	{
		D0 = 15.88;
	}
	if(c == 22 || c == 23 || c == 24 || c == 25 || c == 26 || c == 27 || c == 28 || c == 29 || c == 30 || c == 31)
	{
		D0 = 19.05;
	}
	if(c == 32 || c == 33 || c == 34 || c == 35 || c == 36 || c == 37 || c == 38 || c == 39 || c == 40 || c == 41)
	{
		D0 = 22.23;
	}
	if(c == 42 || c == 43 || c == 44 || c == 45 || c == 46 || c == 47 || c == 48 || c == 49 || c == 50 || c == 51)
	{
		D0 = 25.4;
	}
	if(c == 52 || c == 53 || c == 54 || c == 55 || c == 56 || c == 57 || c == 58 || c == 59 || c == 60 || c == 61)
	{
		D0 = 31.75;
	}
	if(c == 62 || c == 63 || c == 64 || c == 65)
	{
		D0 = 38.1;
	}
	if(c == 66 || c == 67 || c == 68 || c == 69)
	{
		D0 = 50.8;
	}
	if(c == 70 || c == 71 || c == 72)
	{
		D0 = 63.5;
	}
	if(c == 73 || c == 74 || c == 75)
	{
		D0 = 76.2;
	}
	return D0;
}
double CALC_RE(double r, double v, double d, double vi)
{
	return ((r*v*d*pow(10,-3))/vi);
}
double CALC_PR(double vi, double cp, double k)
{
	return ((vi*cp)/k);
}
double CALC_NU(double rea, double pra)
{
	double Nua;
	if(rea <= 2000)
	{
		Nua = 4.36;
	}
	if(rea > 10000 && rea < 5000000 && pra > 0.5 && pra < 2000)
	{
		double f1,f2,f,f3,f4,up,down;
		f1 = (1.58*log(rea) - 3.28);
		f = pow(f1,-2);
		f2 = (f/2);
		f3 = pow(f2,0.5);
		f4 = pow(pra,0.67);
		up = f2*rea*pra;
		down = 1.07 + (12.7*f3*(f4-1));	
		Nua = up/down;
	}
	if(rea > 2300 && rea < 10000 && pra > 0.5 && pra < 2000)
	{
		double f1,f2,f,f3,f4,up,down;
		f1 = (1.58*log(rea) - 3.28);
		f = pow(f1,-2);
		f2 = f/2;
		f3 = pow(f2,0.5);
		f4 = pow(pra,0.67);
		up = f2*(rea-1000)*pra;
		down = 1 + (12.7*f3*(f4-1));
		Nua = up/down;
	}
	return Nua;
}
double CALC_H(double n, double k, double d)
{
	double d1;
	d1=d*pow(10,-3);
	return ((n*k)/d1);
}
double CALC_UC(double d0, double di, double hi, double k, double h0)
{
	double u,u1,f1,f2,f3,f4,f5;
	f1 = (d0/di);
	f2 = (1/hi);
	f3 = (1/h0);
	f4 = (d0/(2*k));
	f5 = log(f1);
	u = ((f1*f2) + (f4*f5) + f3);
	u1 = 1/u;
	return u1; 
}
double CALC_UF(double u, double rf)
{
	double u1,u2,u3;
	u1 = 1/u;
	u2 = u1 + rf;
	u3 = 1/u2;
	return u3;
}
double CALC_CF(double uf, double uc)
{
	double c;
	c = uc/uf;
	return c;
}
double CALC_OD(double c)
{
	double o;
	o = (c - 1)*100;
	return o;
}
double CALC_LMTD(double Tin, double Tout)
{
	double T,t1,t;
	T = (Tin/Tout);
	t1 = Tin - Tout;
	t = t1/(log(T));
	return t;	
}
void Prop_water(double temp)
{
	a1 = fopen(".\\Properties\\water_property.txt","r");
	int i;
	double t[100],rho[100],vis[100],cp[100],kt[100];
	for(i = 0;i<100;i++)
	{
		if(i == 0)
		{
			fscanf(a1,"%*s %*s %*s %*s %*s");
		}
		fscanf(a1,"%lf %lf %lf %lf %lf",&t[i],&rho[i],&vis[i],&cp[i],&kt[i]);
	}
	fclose(a1);
	for(i = 0;i<100;i++)
	{
		if(temp == t[i])
		{
			v = vis[i];
			r = rho[i];
			c = cp[i];
			k = kt[i];
			goto ENTER_VALUE;
		}
		if(temp > t[i-1] && temp < t[i])
		{
			v = vis[i] + (((t[i-1] - temp)*(vis[i+1] - vis[i]))/(t[i+1] - temp)); 
			r = rho[i] + (((t[i-1] - temp)*(rho[i+1] - rho[i]))/(t[i+1] - temp));
			c = cp[i] + (((t[i-1] - temp)*(cp[i+1] - cp[i]))/(t[i+1] - temp)); 
			k = kt[i] + (((t[i-1] - temp)*(kt[i+1] - kt[i]))/(t[i+1] - temp));
			if(temp > 70);
			{
				v = ((((temp - t[i])*(temp - t[i+1]))/((t[i-1] - t[i])*(t[i-1] - t[i+1])))*vis[i-1]) + ((((temp - t[i-1])*(temp - t[i+1]))/((t[i] - t[i-1])*(t[i] - t[i+1])))*vis[i]) + ((((temp - t[i-1])*(temp - t[i]))/((t[i+1] - t[i-1])*(t[i+1] - t[i])))*vis[i+1]);
				r = ((((temp - t[i])*(temp - t[i+1]))/((t[i-1] - t[i])*(t[i-1] - t[i+1])))*rho[i-1]) + ((((temp - t[i-1])*(temp - t[i+1]))/((t[i] - t[i-1])*(t[i] - t[i+1])))*rho[i]) + ((((temp - t[i-1])*(temp - t[i]))/((t[i+1] - t[i-1])*(t[i+1] - t[i])))*rho[i+1]);
				c = ((((temp - t[i])*(temp - t[i+1]))/((t[i-1] - t[i])*(t[i-1] - t[i+1])))*cp[i-1]) + ((((temp - t[i-1])*(temp - t[i+1]))/((t[i] - t[i-1])*(t[i] - t[i+1])))*cp[i]) + ((((temp - t[i-1])*(temp - t[i]))/((t[i+1] - t[i-1])*(t[i+1] - t[i])))*cp[i+1]);
				k = ((((temp - t[i])*(temp - t[i+1]))/((t[i-1] - t[i])*(t[i-1] - t[i+1])))*kt[i-1]) + ((((temp - t[i-1])*(temp - t[i+1]))/((t[i] - t[i-1])*(t[i] - t[i+1])))*kt[i]) + ((((temp - t[i-1])*(temp - t[i]))/((t[i+1] - t[i-1])*(t[i+1] - t[i])))*kt[i+1]);
				goto ENTER_VALUE;
			}
			goto ENTER_VALUE;
		}
	}
	ENTER_VALUE:
	return;
}
void Prop_air(double temp)
{
	a2 = fopen(".\\Properties\\air_property.txt","r");
	int i;
	double t[100],rho[100],vis[100],cp[100],kt[100],p[100];
	for(i = 0;i<100;i++)
	{
		if(i == 0)
		{
			fscanf(a2,"%*s %*s %*s %*s %*s %*s");
		}
		fscanf(a2,"%lf %lf %lf %lf %lf",&t[i],&rho[i],&vis[i],&cp[i],&kt[i],&p[i]);
	}
	fclose(a2);
	for(i = 0;i<100;i++)
	{
		if(temp == t[i])
		{
			v = vis[i];
			r = rho[i];
			c = cp[i];
			k = kt[i];
			goto ENTER_VALUE;
		}
		if(temp > t[i-1] && temp < t[i])
		{
			v = vis[i] + (((t[i-1] - temp)*(vis[i+1] - vis[i]))/(t[i+1] - temp)); 
			r = rho[i] + (((t[i-1] - temp)*(rho[i+1] - rho[i]))/(t[i+1] - temp));
			c = cp[i] + (((t[i-1] - temp)*(cp[i+1] - cp[i]))/(t[i+1] - temp)); 
			k = kt[i] + (((t[i-1] - temp)*(kt[i+1] - kt[i]))/(t[i+1] - temp));
			goto ENTER_VALUE;
		}
	}
	ENTER_VALUE:
	return;
}
void Prop_oil(double temp)
{
	a3 = fopen(".\\Properties\\oil_property.txt","r");
	int i;
	double t[20],rho[20],vis[20],cp[20],kt[20];
	for(i = 0;i<20;i++)
	{
		if(i == 0)
		{
			fscanf(a3,"%*s %*s %*s %*s %*s");
		}
		fscanf(a3,"%lf %lf %lf %lf %lf",&t[i],&vis[i],&rho[i],&cp[i],&kt[i]);
	}
	fclose(a3);
	for(i = 0;i<100;i++)
	{
		if(temp == t[i])
		{
			v = vis[i];
			r = rho[i];
			c = cp[i];
			k = kt[i];
			goto ENTER_VALUE;
		}
		if(temp > t[i-1] && temp < t[i])
		{
			v = ((((temp - t[i])*(temp - t[i+1]))/((t[i-1] - t[i])*(t[i-1] - t[i+1])))*vis[i-1]) + ((((temp - t[i-1])*(temp - t[i+1]))/((t[i] - t[i-1])*(t[i] - t[i+1])))*vis[i]) + ((((temp - t[i-1])*(temp - t[i]))/((t[i+1] - t[i-1])*(t[i+1] - t[i])))*vis[i+1]);
			r = ((((temp - t[i])*(temp - t[i+1]))/((t[i-1] - t[i])*(t[i-1] - t[i+1])))*rho[i-1]) + ((((temp - t[i-1])*(temp - t[i+1]))/((t[i] - t[i-1])*(t[i] - t[i+1])))*rho[i]) + ((((temp - t[i-1])*(temp - t[i]))/((t[i+1] - t[i-1])*(t[i+1] - t[i])))*rho[i+1]);
			c = ((((temp - t[i])*(temp - t[i+1]))/((t[i-1] - t[i])*(t[i-1] - t[i+1])))*cp[i-1]) + ((((temp - t[i-1])*(temp - t[i+1]))/((t[i] - t[i-1])*(t[i] - t[i+1])))*cp[i]) + ((((temp - t[i-1])*(temp - t[i]))/((t[i+1] - t[i-1])*(t[i+1] - t[i])))*cp[i+1]);
			k = ((((temp - t[i])*(temp - t[i+1]))/((t[i-1] - t[i])*(t[i-1] - t[i+1])))*kt[i-1]) + ((((temp - t[i-1])*(temp - t[i+1]))/((t[i] - t[i-1])*(t[i] - t[i+1])))*kt[i]) + ((((temp - t[i-1])*(temp - t[i]))/((t[i+1] - t[i-1])*(t[i+1] - t[i])))*kt[i+1]);
			goto ENTER_VALUE;
		}
	}
	ENTER_VALUE:
	return;
}
double calc_f(double re, double K, double D)
{
	double K1,A,B,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,C;
	K1 = K/D;
	u1 = (K1/3.827);
	u2 = (K1/7.7918);
	u3 = pow(u2,0.9924);
	u4 = 5.3326/(208.815 + re);
	u5 = pow(u4,0.9345);
	u6 = u3 + u5;
	u7 = log(u6);
	u8 = 4.567/re;
	u9 = u8*u7;
	u10 = u1 - u9;
	A = log(u10);
	u11 = 5.0272/re;
	u12 = K1/3.7065;
	u13 = u11*A;
	B = u12 - u13;
	C = log(B);
	u14 = -2*C;
	u15 = 1/u14;
	u16 = pow(u15,2);
	return u16;
}
double calc_pressuredrop_Tubein(double f, double n, double l, double d1, double u, double rho)
{
	double a,b,c,d,e;
	a = 2*l/d1;
	b = n*50;
	c = pow(u,2);
	d = a+b;
	e = 4*f*d*0.5*rho*c;
	return e;
	
}
double calc_pressuredrop_Annulus(double f, double n, double l, double d1, double u, double rho,double h)
{
	double a,b,c,d,e,g,z;
	a = 2*l/d1;
	b = n*75;
	c = pow(u,2);
	z = calc_90(h);
	g = (2*(n-1)*z);
	d = a+b+g;
	e = 4*f*d*0.5*rho*c;
	return e;
}
double calc_90(double a)
{
	double l,up,down,mul,ans=0;
	double c[5],e[5];
	int i,j;
	a6 = fopen(".\\PressureDrop\\bend.txt","r");
	for(i = 0;i < 5;i++)
	{
		fscanf(a6,"%*s%lf \t %*s%lf",&c[i],&e[i]);
	}
	fclose(a6);
	for(i=0;i<5;i++)
	{
		mul=1;
		for(j=0;j<5;j++)
		{
			up=1;
			down=1;
			if(j!=i)
			{
				up=up*(a-c[j]);
				down=down*(c[i]-c[j]);
			}
			mul=mul*(up/down);
		}
		ans = ans + (mul * e[i]);
	}
	l = ans;
	for (i = 0; i< 5;i++)
	{
		if(a == c[i])
		{
			l = e[i];
		}
	}
	return l;
}
void makefile()
{
	fp = fopen(".\\Input\\Flow_input.txt","r");
	fp1 = fopen(".\\Output\\Hotside_Calculation.txt","w");
	fp2 = fopen(".\\Output\\Coldside_calculation.txt","w");
	fp3 = fopen(".\\Output\\Tubeside_Calculation.txt","w");
	fp4 = fopen(".\\Output\\Annulusside_Calculation.txt","w");
	fp5 = fopen(".\\Selection Option\\Tube_side_selection_Option.txt","w");
	fp6 = fopen(".\\Selection Option\\Annulus_tube_selection_option.txt","w");
	fp7 = fopen(".\\Output\\General_calculation.txt","w");
	fp8 = fopen(".\\Output\\Epsilon_Variation_output.xlsx","w");
	fp9 = fopen(".\\Output\\Pressure_Drop_Calculation_different_length.xlsx","w");
	return;
}
