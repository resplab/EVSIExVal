#include <Rcpp.h>
using namespace Rcpp;

#define max(a,b)            \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b);  \
  _a > _b ? _a : _b; })     \


#define min(a,b)              \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b);    \
  _a > _b ? _b : _a; })       \



//Principles: all preprocessing and checks etc are done by the R wrapper. Expects a numeric matrix with first three columsn being prev, se, sp.
// [[Rcpp::export]]
List CEVSI(NumericMatrix samples, double z, NumericVector futureSampleSizes, int nSim, bool debug=false)
{
  int nSamples=samples.nrow();
  int nNStars=futureSampleSizes.length();
  int nStars[nNStars]; //Will copy futureSampleSizes for faster operations
  double data[nSamples][13]; //(double *)malloc(nSamples*(3+6+3+1)*sizeof(double)); //3 original draws, 6 log(theta) and log(1-theta), 3 NBs, and 1 W
  long double lw[nSamples];
  double ENBs[3]={0,0,0};
  double EVPI=0;
  double EVSI[nNStars];

  List out=List::create();

  for(int i=0;i<nSamples;i++)
  {
    data[i][0]=samples(i,0); //prev
    data[i][1]=samples(i,1); //sp
    data[i][2]=samples(i,2); //sp

    data[i][3]=std::log(data[i][0]);
    data[i][4]=std::log(1-data[i][0]);
    data[i][5]=std::log(data[i][1]);
    data[i][6]=std::log(1-data[i][1]);
    data[i][7]=std::log(data[i][2]);
    data[i][8]=std::log(1-data[i][2]);

    data[i][9]=0;
    data[i][10]=data[i][0]*data[i][1]-(1-data[i][0])*(1-data[i][2])*z/(1-z);
    data[i][11]=data[i][0]-(1-data[i][0])*z/(1-z);

    data[i][12]=0;  //Will be the log(weight)
  }

  for(int i=0;i<nNStars;i++)
  {
    nStars[i]=futureSampleSizes(i);
  }

  for(int iStar=0;iStar<nNStars;iStar++)  EVSI[iStar]=0;

  for(int iSample=0;iSample<nSamples;iSample++) //this is the truth!
  {
    double trueNBs[]={data[iSample][9], data[iSample][10], data[iSample][11]};
    double truePrev=data[iSample][0];
    double trueSe=data[iSample][1];
    double trueSp=data[iSample][2];

    EVPI+=max(trueNBs[0],max(trueNBs[1],trueNBs[2]));
    ENBs[0]+=trueNBs[0];
    ENBs[1]+=trueNBs[1];
    ENBs[2]+=trueNBs[2];

    for(int iSim=0;iSim<nSim;iSim++)
    {
      for(int iStar=0;iStar<nNStars;iStar++)
      {
        int n=nStars[iStar];

        //1. Generate binomial numbers
        int nD=R::rbinom(n,truePrev);
        int nTP=R::rbinom(nD,trueSe);
        int nTN=R::rbinom(n-nD,trueSp);

        //2. Weigh them against sample and calculate log weights
        long double minLw, maxLw;
        for(int jSample=0;jSample<nSamples;jSample++)
        {
          long double _lw=(long double)nD*data[jSample][3]+(n-nD)*data[jSample][4]+
            (long double)nTP*data[jSample][5]+(nD-nTP)*data[jSample][6]+
            (long double)nTN*data[jSample][7]+(n-nD-nTN)*data[jSample][8];
          if(jSample==0)
          {
            minLw=_lw; maxLw=_lw;
          }
          else
          {
            minLw=min(_lw,minLw); maxLw=max(_lw,maxLw);
          }
          lw[jSample]=_lw;
        }//jSample

        //3. Back transform weights and do the weighted averaging of NBs
        long double NBmodel=0, NBall=0;
        for(int jSample=0;jSample<nSamples;jSample++)
        {
          long double _w=std::exp2l(lw[jSample]);
          //Rprintf("%f\n",w);
          NBmodel=NBmodel+_w*(long double)data[jSample][10];
          NBall=NBall+_w*(long double)data[jSample][11];
        }
        int winner=0;
        if(NBmodel>0)
          if(NBall>NBmodel) winner=2; else winner=1;
        else
          if(NBall>0) winner=2;

        EVSI[iStar]+=trueNBs[winner];
      } //iStar
    } //iSim
  }//iSample

  double EVCI=max(ENBs[0],max(ENBs[1],ENBs[2]))/nSamples;
  EVPI=EVPI/nSamples-EVCI;
  for(int iStar=0;iStar<nNStars;iStar++)  EVSI[iStar]=EVSI[iStar]/(nSim*nSamples)-EVCI;

  if(debug)
  {
    out["nSamples"]=nSamples;
    out["nNStars"]=nNStars;
    out["data"]=transpose(NumericMatrix(13, nSamples, &data[0][0]));
  }

  out["EVPI"]=EVPI;
  NumericVector tmp(nNStars);
  for(int i=0;i<nNStars;i++) tmp(i)=EVSI[i];
  out["EVSI"]=tmp;

  return out;
}




