Curive fitting routing
originally by Jørgen, adopted by Fabio

Python script to fit 2 Generalized Lorentzian (GLO) or (EGLO) functions to (gamma,xn) data, some pygmy resonances of Standard Lorentzian form, and the scissorrs mode also of Standard Lorentzian form.

# Remarks:
The fitting routine doues not seem to able to determine the temperature. This makes sense, as we have too little data (in my example at least!) to determine the behaviour at very low energies. Note however, that the fitted values depend slightly on the chosen temp. It can be chosen by the hack of setting the limits Tmin and Tmax in the script very close to the good-guess temperature.
Also, we noticed that the standard deviation e.g. for the GLO witdth is quite large. However ,lLooking at the relativly poor quality of the data data vs the very high number of free parameters I think that makes a lot of sense.

# Dependencies:
* matplotlib.pyplot
* python-tk

The repo include the script itself and some example data. To run it, simply type

`python curve_fit.py`


#Output:
```
Shall the (middel/otherwise defined) part of the OCL_data be skipped? (y/n)
n
choice: take all data
Fit results:         Starting value      Optimal value   Standard deviation
Parameter     (E)GLO1_E:             11.30             11.37             0.43
Parameter (E)GLO1_gamma:             3.20              3.75              1.39
Parameter (E)GLO1_sigma:             290.00            310.10            50.47
Parameter     (E)GLO2_E:             14.15             14.35             0.21
Parameter (E)GLO2_gamma:             5.50              4.73              1.19
Parameter (E)GLO2_sigma:             340.00            303.80            57.37
Parameter      (E)GLO_T:             0.34              0.30              0.16
Parameter        SLO1_E:             1.50              1.52              0.08
Parameter    SLO1_gamma:             0.50              0.28              0.55
Parameter    SLO1_sigma:             0.68              0.08              0.06
Parameter        SLO2_E:             2.50              2.13              0.06
Parameter    SLO2_gamma:             0.50              0.81              0.28
Parameter    SLO2_sigma:             0.35              0.70              0.15
Parameter        SLO3_E:             7.50              7.77              2.72
Parameter    SLO3_gamma:             5.45              6.04              4.07
Parameter    SLO3_sigma:             20.00             20.08             14.52

 (E)GLO temp:    Start   Tmin    Tmax    Opt     Std
         0.34    0.30    0.40    0.30 +- 0.16
```
![Fit Results](https://github.com/fzeiser/curve_fit_gamma_strength/blob/master/fit.png "Fit Results")


# Preparation:
In order to get the input files with the cross sections converted to gamma-ray strength, you might run something like this in root
```cpp

const float factor = 8.674E-08; // const. factor in mb^(-1) MeV^(-2)
int    i, j;

ifstream strengthfile("../../strength.nrm");
int nLines=0;
//std::ifstream myfile("main.cpp");
std::string lines;
while (std::getline(strengthfile, lines))
       ++nLines;
int nPoints = nLines/2;
//cout << nPoints << endl;
ofstream outfile("AA_gSF_240Pu_strength.dat");
outfile << "Eg(MeV)" << "\t" <<  "f(MeV^-3)" << "\t" << "f_Err(MeV^-3)" << endl;

ifstream strengthfile("../../strength.nrm");
    float strength[90],strengtherr[90],energy[90],energyerr[90];
    int i = 0;
   float a0 =  -0.8360; // adopt to your data!
   float a1 =   0.1280; // adopt to your data!
    float x;    
    while(strengthfile){
        strengthfile >> x;

        //"initialize"
        strength[i]=0;
        energy[i]=0;
        energyerr[i]=0;
        strengtherr[i]=0;

        if(i<nPoints){
            strength[i] = x;
            energy[i] = a0 + (a1*i);
            energyerr[i] = 0.0;
           cout << i << "\t" << energy[i] << endl;
        }   
        else{strengtherr[i-nPoints] = x;
             cout << i << "\t" << strengtherr[i-nPoints] << endl; }

        i++;
    }
    
    for(i = 0; i < nPoints; i++){
            if ( strength[i]>0 && energy[i]>0){ 
        outfile << energy[i] << "\t" << strength[i]  << "\t" << strengtherr[i] << endl;
            }
        }

    outfile.close();
    cout << "\n strengthfile was written \n" << endl;
    TGraphErrors *strengthexp = new TGraphErrors(nPoints,energy,strength,energyerr,strengtherr);
```

and for the cross-sections

```cpp

    ///////////////////////////////////////////////////////
    
    float gdr[70],gdr_err[70],gdr_energy[70],gdr_energyerr[70];
    i = 0;
    ifstream gdrfile("239Pu_gurevich_1976_g_abs.txt");
    
    ofstream outfile("AA_gSF_239Pu_gurevich_1976_g_abs.txt");
    outfile << "Eg(MeV)" << "\t" <<  "f(MeV^-3)" << "\t" << "f_Err(MeV^-3)" << endl;
    
    while(gdrfile){
        gdrfile >> x >> y >> z;
//      gdrfile >> x >> y;            // errors not read in because non existen
        if(i<47){
            gdr_energy[i]   = x;
            gdr[i]          = y;
            gdr_err[i]      = z;
//          gdr_err[i]      = 0; // errors not read in because non existen
            gdr_energyerr[i]= 0;
        }
        i++;
    }
    for(i = 0; i < 47; i++){
        gdr[i]     = factor*gdr[i]    /gdr_energy[i];
        gdr_err[i] = factor*gdr_err[i]/gdr_energy[i];

    if ( gdr_energy[i]>0 && gdr[i]>0){  
        outfile << gdr_energy[i] << "\t" << gdr[i]  << "\t" << gdr_err[i] << endl;
        }
    }
//    gdr[0]=0; gdr[1]=0; gdr[2]=0;
//    gdr_energy[0]=20; gdr_energy[1]=20; gdr_energy[2]=20;
    outfile.close();
    cout << "\n strengthfile was written \n" << endl;

    TGraphErrors *gdrexp1 = new TGraphErrors(48,gdr_energy,gdr,gdr_energyerr,gdr_err);

    ///////////////////////////////////////////////////////
```
