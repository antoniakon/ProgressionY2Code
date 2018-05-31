import java.io.File

/**
  * Created by Antonia Kontaratou.
  * Saturated model for main effects as random effects by putting a gamma prior on the precision of the effects.
  * This version uses the new notation with alpha and beta to be consistent with the maths.
  */
object MainEffRandom {

  //---------------------- Main effects ---------------------------//
  // the model is: Xijk|mu,aj,bk, tau~N(mu+aj+bk,tau^-1)

  def saturatedModelCatMainEffs(noOfIter: Int, thin:Int, y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], nj: Int, nk: Int, structure: DVStructure, alphaPriorMean: Double, betaPriorMean: Double, mu0: Double, tau0: Double, a: Double, b:Double, aPrior: Double, bPrior:Double) = {

    val N = y.length //the number of observations
    val sampleNo= noOfIter/thin+1 //the number of samples created from the MCMC
    val mat_mt = new DenseMatrix[Double](sampleNo, 2) //to store mu and tau
    val alphaCoefs = new DenseMatrix[Double](sampleNo, nj) //to store the coefficients for variable alpha
    val betaCoefs = new DenseMatrix[Double](sampleNo, nk) // to store the coefficients for variable beta

    val curAlpha = DenseVector.zeros[Double](nj) //the current values of the coefficients for alpha at every iteration. Initialised with 0.
    val curBeta = DenseVector.zeros[Double](nk) //the current values of the coefficients for alpha at every iteration. Initialised with 0.

    var mu = 0.0 //initialise the sampler at mu=0
    var tau = 1.0 // initialise the sampler at tau=1.0

    //Initialise with 0 in the first row the matrices where we will store the values to start the Gibbs sampler from 0s
    mat_mt(0,::):= DenseVector(mu,tau).t
    alphaCoefs(0,::):= curAlpha.t
    betaCoefs(0,::):= curBeta.t
    val SumX= y.toArray.sum // the sum of the values of all the observations

    var tauAlpha= 1.0 //initialise the precision for alpha (this is actually tau_alpha)
    var tauBeta= 1.0 //initialise the precision for beta (this is actually tau_beta)

    var sumaj=0.0 //this will be used for sampling the precision of alpha from the FCD. Here it is initialised at 0.
    var sumbk=0.0 //this will be used for sampling the precision of beta from the FCD. Here it is initialised at 0.

    var ind= 0 //index used for thinning

    for(i<-1 until noOfIter ){
      println("iter :" + i)
      for (j<-0 until nj){
        sumaj= sumaj+pow((curAlpha(j)-alphaPriorMean),2) //estimate the sum used in sampling from Gamma distribution for the precision of alpha
      }

      for (k<-0 until nk){
        sumbk= sumbk+pow((curBeta(k)-betaPriorMean),2) //estimate the sum used in sampling from Gamma distribution for the precision of beta
      }

      tauAlpha= breeze.stats.distributions.Gamma(aPrior+nj/2.0, 1.0/(bPrior+0.5*sumaj)).draw() //sample the precision of alpha from gamma
      tauBeta= breeze.stats.distributions.Gamma(aPrior+nk/2.0, 1.0/(bPrior+0.5*sumbk)).draw() // sample the precision of beta from gamma

      //Update mu and tau
      //mu = 3.10
      //tau = 1.0
      val varMu= 1.0/(tau0+N*tau) //the variance for mu
      val meanMu= (mu0*tau0+tau*(SumX - sumAllEffects(y,alpha, beta, curAlpha, curBeta, nj, nk)))*varMu
      mu= breeze.stats.distributions.Gaussian(meanMu,sqrt(varMu)).draw()
      tau= breeze.stats.distributions.Gamma(a+N/2.0, 1.0/(b+0.5*YminusMuAndEffects(y,alpha,beta,mu,curAlpha,curBeta))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β

      // Update alphaj
      for(j <-0 until nj){
        val SXalphaj= structure.calcAlphaSum(j)  // the sum of the observations that have alpha==j
        val Nj= structure.calcAlphaLength(j)     // the number of the observations that have alpha==j
        val SumBeta= sumBetaEffGivenAlpha(y, alpha, beta, j, curBeta)
        val varPalpha = 1.0/(tauAlpha + tau*Nj)   //the variance for alphaj
        val meanPalpha = (alphaPriorMean*tauAlpha+ tau*(SXalphaj-Nj*mu-SumBeta))*varPalpha
        curAlpha(j) = breeze.stats.distributions.Gaussian(meanPalpha, sqrt(varPalpha)).draw()
      }
      //println("iter: "+ i+ curAlpha)

      //Update betak
      for(k <-0 until nk){
        val SXbetak= structure.calcBetaSum(k)   // the sum of the observations that have beta==k
        val Nk= structure.calcBetaLength(k)   // the number of the observations that have beta==k
        val SumAlpha= sumAlphaEffGivenBeta(y,alpha,beta, k, curAlpha)
        val varPbeta = 1.0/(tauBeta+ tau*Nk)
        val meanPbeta = (betaPriorMean*tauBeta + tau*(SXbetak-Nk*mu-SumAlpha))*varPbeta
        curBeta(k)= breeze.stats.distributions.Gaussian(meanPbeta, sqrt(varPbeta)).draw()
        }
      //println("iter: "+ i+ curBeta)
      // Before using thinning.
      //mat_mt(i,::):= DenseVector(mu,tau).t
      //alphaCoefs(i,::):= curAlpha.t
      //betaCoefs(i,::):= curBeta.t

      // After using thinning
      //curAlpha:= DenseVector(-1.0,4.0,-3.0)
      //curBeta:= DenseVector(-1.0,3.5,-2.0,-0.5)
      if((i%thin).equals(0)){
        //println("in if")
        ind=ind+1
        mat_mt(ind, ::) := DenseVector(mu, tau).t
        alphaCoefs(ind, ::) := curAlpha.t
        betaCoefs(ind, ::) := curBeta.t
      }
      sumaj=0.0
      sumbk=0.0

    }
    (mat_mt,alphaCoefs,betaCoefs)

  }
  /**
    * Add all the beta effects for a given alpha.
    */
  def sumBetaEffGivenAlpha(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaIndex: Int, betaEff: DenseVector[Double]):Double={
    val N= y.length
    var sum=0.0
    for(i <- 0 until N){                    // run through all the observations
      if(alpha(i)== alphaIndex+1){          // if alpha of the current observation == current alphaIndex [+1 because of the difference in the dataset notation (alpha=1,2,...) and Scala indexing that starts from 0]
        sum= sum + betaEff(beta(i)-1)       // add to the sum the current effect of the observation's beta (-1 for consistency with the Scala vector indexing again)
      }
    }
    sum
  }
  /**
    * Add all the alpha effects for a given beta.
    */
  def sumAlphaEffGivenBeta(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], betaIndex: Int, alphaEff: DenseVector[Double]):Double={
    val N= y.length
    var sum=0.0
    for(i <- 0 until N){                    // run through all the observations
      if(beta(i)== betaIndex+1){            // if beta of the current observation == current betaIndex [+1 because of the difference in the dataset notation (beta=1,2,...) and Scala indexing that starts from 0]
        sum= sum + alphaEff(alpha(i)-1)     // add to the sum the current effect of the observation's alpha (-1 for consistency with the Scala vector indexing again)
      }
    }
    sum
  }

  /**
    * Calculate the Yi-mu-alpha_eff-beta_eff. To be used in estimating tau
    */
  def YminusMuAndEffects(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], mu: Double, alphaEff: DenseVector[Double], betaEff: DenseVector[Double]):Double={
    val N= y.length
    val YminusMuNatUniEff = DenseVector.zeros[Double](N)
    for (i <- 0 until N){
      YminusMuNatUniEff(i)= y(i)-alphaEff(alpha(i)-1)-betaEff(beta(i)-1)-mu //run through all the observations and from the response subtract (mu + the effect of its alpha + the effect of its beta)
    }
    sqr(YminusMuNatUniEff).toArray.sum //square all the elements of the vector and find the sum
  }

  /**
    * Calculate the sum of all the alpha and all the beta effects for all the observations.
    */
  def sumAllEffects(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaEff: DenseVector[Double], betaEff: DenseVector[Double], nj: Int, nk: Int):Double={
    var sumAlpha=0.0
    var sumBeta=0.0
    // For the effects of alpha
    for(i <- 0 until nk){   // run though all the betas to find the alphas' effects.
      sumAlpha= sumAlpha + sumAlphaEffGivenBeta(y,alpha,beta, i, alphaEff)
    }
    // For the effects of beta
    for(i <- 0 until nj){   // run though all the alphas to find the betas' effects.
      sumBeta= sumBeta + sumBetaEffGivenAlpha(y,alpha,beta, i, betaEff)
    }
    sumAlpha+sumBeta  //their sum is the sum of all the alphas' + betas' effects
  }

  def sqr(x: DenseVector[Double]) = x * x

  def time[A](f: => A) = {
    val s = System.nanoTime
    val ret = f
    println( "time: " +(System.nanoTime-s)/1e6+ "ms")
    ret
  }

  def main(args: Array[String]): Unit = {
    // Read the data
    val data = csvread(new File("./simulNoInter.csv"))
    val sampleSize = data.rows
    val y = data(::, 0)
    val alpha = data(::, 1).map(_.toInt)
    val beta = data(::, 2).map(_.toInt)
    val structure = new DVStructure(y, alpha, beta)

    // Parameters
    val nj = data(::, 1).toArray.distinct.length
    val nk = data(::, 2).toArray.distinct.length
    val noOfIters = 1000000
    val thin = 100
    val aPrior = 1
    val bPrior = 0.0001
    val alphaPriorMean = 0.0
    val alphaPriorTau = 0.0001
    val betaPriorMean = 0.0
    val betaPriorTau = 0.0001
    val mu0 = 0.0
    val tau0 = 0.0001
    val a = 1
    val b = 0.0001

    /**
      * SATURATED MODEL. MAIN EFFECTS AS RANDOM.
      * Saturated model main effects with random effects and gamma prior on the precision of the effects.
      */

    val (test_mt,alpha_est,beta_est)= time(saturatedModelCatMainEffs(noOfIters, thin, y, alpha, beta, nj, nk, structure, alphaPriorMean, betaPriorMean, mu0, tau0, a,b, aPrior, bPrior))
    val mt= mean(test_mt(::,*)).t
    val alphaEstim= mean(alpha_est(::,*)).t
    val betaEstim= mean(beta_est(::,*)).t

    println("Results: ")
    println("mu and tau:" + mt)
    println("Estimates for alpha: "+ alphaEstim)
    println("Estimates for beta: "+ betaEstim)

    // The expected values are:
    // mu = 3.1, tau = 1.0
    // For alpha:= DenseVector(-1.0,4.0,-3.0)
    // For beta:= DenseVector(-1.0,3.5,-2.0,-0.5)

    //Plot Results
    val resMat1= DenseMatrix.horzcat(test_mt,alpha_est, beta_est)
    val plot1= new PlotResults
    plot1.plotResults(resMat1, List("Mean","tau", "alpha 1", "alpha 2", "alpha 3","beta 1","beta 2","beta 3","beta 4"))

  }
}
