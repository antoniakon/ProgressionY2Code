import java.io.File

/**
  * Created by Antonia Kontaratou.
  * Saturated model for Main + interaction effects as fixed effects using sum to 0 constraint to solve the identifiability problem.
  */
object MainAndInterEffSumToZero {
  def saturatedModelCatMainEfAndInteractions(noOfIter: Int, y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], lalpha: Int, lbeta: Int, structure: DVStructure, ujPriorMean: Double, ujPriorVar: Double, nkPriorMean: Double, nkPriorVar: Double, muPrior: Double, varPrior: Double, a: Double, b:Double, interPriorMean: DenseMatrix[Double], interPriorCov: DenseMatrix[Double]) = {

    val N = y.length
    val mat_mt = new DenseMatrix[Double](noOfIter, 2) //to store mu and tau
    val alphaCoefs = new DenseMatrix[Double](noOfIter, lalpha) //to store the coefficients for alpha
    val betaCoefs = new DenseMatrix[Double](noOfIter, lbeta) // to store the coefficients for beta
    val interCoefs = new DenseMatrix[Double](noOfIter, lbeta * lalpha) // to store the coefficients for the interactions

    val curAlpha = DenseVector.zeros[Double](lalpha) //the current values of the coefficients for every iteration. Initialised with 0.
    val curBeta = DenseVector.zeros[Double](lbeta)
    val curInter = DenseMatrix.zeros[Double](lalpha, lbeta)
    val SYi = y.toArray.sum

    var mu = 0.0 //initialise the sampler at a=0
    var tau = 1.0 // initialise the sampler at tau=0

    //Initialise with 0 in the first row the matrices where we will store the values to start the Gibbs sampler from 0s
    mat_mt(0, ::) := DenseVector(mu, tau).t
    alphaCoefs(0, ::) := curAlpha.t
    betaCoefs(0, ::) := curBeta.t
    interCoefs(0, ::) := curInter.t.toDenseVector.t

    for (i <- 1 until noOfIter) {

      //Update mu and tau
      //mu = 3.10
      tau = 1.0
      val varMu = 1.0 / ((1.0 / varPrior) + N * tau)
      val meanMu = ((muPrior / varPrior) + tau * (SYi - sumAllMainInterEff(y, alpha, beta, curAlpha, curBeta, lalpha, lbeta, curInter))) / (1.0 / varMu)
      mu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
      tau = breeze.stats.distributions.Gamma(a + N / 2.0, 1.0 / (b + YminusMainInterEffects(y, alpha, beta, mu, curAlpha, curBeta, curInter) / 2.0)).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β

      //      // Update alphaj
      for (jj <- 0 until lalpha-1) {
        val SYiuj = structure.calcAlphaSum(jj)
        val Nj = structure.calcAlphaLength(jj)
        val Sn = sumBetaEff(y, alpha, beta, jj, curBeta)
        val SinterAlpha = sumInterEffAlpha(y, alpha, beta, jj, curInter)
        val varPU = (1.0 / ((1.0 / ujPriorVar) + tau * Nj))
        val meanPU = ((ujPriorMean / ujPriorVar) + tau * (SYiuj - Nj * mu - Sn - SinterAlpha)) / (1.0 / varPU)
        curAlpha(jj) = breeze.stats.distributions.Gaussian(meanPU, sqrt(varPU)).draw()
      }
      curAlpha(lalpha-1)= -sum(curAlpha(0 until lalpha-1))
      //      // Update betak
      for (jk <- 0 until lbeta-1) {
        val SYink = structure.calcBetaSum(jk)
        val Nk = structure.calcBetaLength(jk)
        val Su = sumAlphaEff(y, alpha, beta, jk, curAlpha)
        val SinterBeta = sumInterEffBeta(y, alpha, beta, jk, curInter)
        val varPN = (1.0 / ((1.0 / nkPriorVar) + tau * Nk))
        val meanPN = ((nkPriorMean / nkPriorVar) + tau * (SYink - Nk * mu - Su - SinterBeta)) / (1.0 / varPN)
        curBeta(jk) = breeze.stats.distributions.Gaussian(meanPN, sqrt(varPN)).draw()
      }
      curBeta(lbeta-1)= -sum(curBeta(0 until lbeta-1))
      // Update Interaction terms
      for (jj <- 0 until lalpha-1) {
        for (jk <- 0 until lbeta-1) {
          val Njk= structure.getDVList(jj,jk).length
          val SYijk = structure.getDVList(jj,jk).sum
          val varPI = 1.0 / ((1.0 / interPriorCov(jj,jk)) + tau * Njk)
          val meanPI=((interPriorMean(jj,jk) / interPriorCov(jj,jk)) + tau * (SYijk - Njk * mu - Njk*(curAlpha(jj)+curBeta(jk)))) / (1.0 / varPI)
          curInter(jj,jk) = breeze.stats.distributions.Gaussian(meanPI, sqrt(varPI)).draw()
        }
      }

      for ( j <-0 until lbeta){
        curInter(lalpha-1,j)= -sum(curInter(0 until lalpha-1,j))
      }

      for ( i <-0 until lalpha){
        curInter(i, lbeta-1)= -sum(curInter(i,0 until lbeta-1))
      }
      //      curAlpha:= DenseVector(-1.0,4.0,-3.0)
      //      curBeta:= DenseVector(-1.0,3.5,-2.0,-0.5)
      //curInter:= DenseMatrix((-0.5,-2.0,0.0,2.5),  (0.0,-0.5,-1.0,1.5), (0.5,2.5,1.0,-4.0))
      mat_mt(i, ::) := DenseVector(mu, tau).t
      alphaCoefs(i, ::) := curAlpha.t
      betaCoefs(i, ::) := curBeta.t
      interCoefs(i, ::) := curInter.t.toDenseVector.t
    }
    (mat_mt, alphaCoefs, betaCoefs, interCoefs)
  }

  /**
    * Add all the beta effects for a given alpha.
    */
  def sumBetaEff(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaIndex: Int, betaEff: DenseVector[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if(alpha(i)== alphaIndex+1){
        sum= sum + betaEff(beta(i)-1)
      }
    }
    sum
  }
  /**
    * Add all the alpha effects for a given beta.
    */
  def sumAlphaEff(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], betaIndex: Int, alphaEff: DenseVector[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if(beta(i)== betaIndex+1){
        sum= sum + alphaEff((alpha(i)-1))
      }
    }
    sum
  }
  /**
    * Calculate the sum of all the alpha and all the beta effects for all the observations.
    */
  def sumAllMainInterEff(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaEff: DenseVector[Double], betaEff: DenseVector[Double], lalpha: Int, lbeta: Int, interEff: DenseMatrix[Double]):Double={
    var sumBeta=0.0
    var sumAlpha=0.0
    var sumInter=0.0
    // For alpha effects
    for(i <- 0 until lbeta){
      sumAlpha= sumAlpha + sumAlphaEff(y,alpha,beta, i, alphaEff)
    }
    sum
    // For beta effects
    for(i <- 0 until lalpha){
      sumBeta= sumBeta + sumBetaEff(y,alpha,beta, i, betaEff)
    }
    // For Interaction effects
    for(i <- 0 until lalpha){
      for(j <- 0 until lbeta){
        sumInter= sumInter + sumInterEff(y,alpha,beta, i, j, interEff)
      }
    }

    sumAlpha+sumBeta+sumInter
  }

  /**
    * Add all the interaction effects for a given beta and alpha.
    */
  def sumInterEff(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaIndex: Int, betaIndex: Int, interEff: DenseMatrix[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if((alpha(i)== alphaIndex+1)&(beta(i)== betaIndex+1)){
        sum = sum + interEff((alpha(i)-1), (beta(i)-1))
      }
    }
    sum
  }

  /**
    * Add all the interaction effects for a given alpha.
    */
  def sumInterEffAlpha(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaIndex: Int, interEff: DenseMatrix[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if((alpha(i)== alphaIndex+1)){
        sum = sum + interEff((alpha(i)-1), (beta(i)-1))
      }
    }
    sum
  }
  /**
    * Add all the interaction effects for a given beta.
    */
  def sumInterEffBeta(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], betaIndex: Int, interEff: DenseMatrix[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if((beta(i)== betaIndex+1)){
        sum = sum + interEff((alpha(i)-1), (beta(i)-1))
      }
    }
    sum
  }
  /**
    * Calculate the Yi-mu-u_eff-n_eff- inter_effe. To be used in estimating tau
    */
  def YminusMainInterEffects(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], mu: Double, alphaEff: DenseVector[Double],betaEff: DenseVector[Double], interEff: DenseMatrix[Double]):Double={
    val N= y.length
    val YminusMuBetaAlphaEff = DenseVector.zeros[Double](N)
    for (i <- 0 until N){
      YminusMuBetaAlphaEff(i)= y(i)-mu-alphaEff(alpha(i)-1)-betaEff(beta(i)-1)-interEff(alpha(i)-1, beta(i)-1)
    }
    sqr(YminusMuBetaAlphaEff).toArray.sum
  }

  /**
    * Add all the beta effects for a given alpha.
    */
  def sumEffAlphaBeta(structure: DVStructure, alphaEff: DenseVector[Double], betaEff: DenseVector[Double], alphaIndex: Int, betaIndex: Int):Double={
    structure.getDVList(alphaIndex, betaIndex).length

  }
  def sqr(x: DenseVector[Double]) = x * x

  def main(args: Array[String]): Unit = {
    // Read the data
    val data = csvread(new File("./simulInter.csv"))

    val sampleSize = data.rows
    val y = data(::, 0)
    val alpha = data(::, 1).map(_.toInt)
    val beta = data(::, 2).map(_.toInt)
    val structure = new DVStructure(y, alpha, beta)

    // Parameters
    val lalpha = data(::, 1).toArray.distinct.length
    val lbeta = data(::, 2).toArray.distinct.length
    val noOfIters = 100000
    val aPrior = 1
    val bPrior = 0.0001
    val ujPriorMean = 0.0
    val ujPriorVar = 1 / 0.0001
    val nkPriorMean = 0.0
    val nkPriorVar = 1 / 0.0001
    val muPrior = 0.0
    val varPrior = 1.0 / 0.0001
    val a = 1
    val b = 0.0001

    val interPriorMeanS = DenseMatrix.fill(lalpha, lbeta){0.0}
    val interPriorCovS= DenseMatrix.fill(lalpha, lbeta){1.0/0.0001} // For random effects sample from Gaussian here

    val (test_mtInterS0,alpha_estInterS0,beta_estInterS0, inter_estS0)= saturatedModelCatMainEfAndInteractions(noOfIters, y, alpha, beta, lalpha, lbeta, structure, ujPriorMean, ujPriorVar, nkPriorMean, nkPriorVar, muPrior, varPrior, a, b, interPriorMeanS, interPriorCovS)
    val mtInterS0= mean(test_mtInterS0(::,*)).t
    val alphaEstimInterS0= mean(alpha_estInterS0(::,*)).t
    val betaEstimInterS0= mean(beta_estInterS0(::,*)).t
    val interEstimInterS0= mean(inter_estS0(::,*)).t

    println("Results: ")
    println("mu and tau:" + mtInterS0)
    println("Estimates for alpha: "+ alphaEstimInterS0)
    println("Estimates for beta "+ betaEstimInterS0)
    println("Estimates for interactions: "+ interEstimInterS0)

    //Plot Results
    val resMat= DenseMatrix.horzcat(test_mtInterS0,alpha_estInterS0, beta_estInterS0, inter_estS0)
    val plot= new PlotResults
    plot.plotResults(resMat, List("Mean", "tau", "alpha 1", "alpha 2", "alpha 3", "beta 1", "beta 2", "beta 3", "beta 4", "a1-b1", "a1-b2", "a1-b3", "a4-b4", "a2-b1", "a2-b2", "a2-b3", "a2-b4", "a3-b1", "a3-b2", "a3-b3", "a3-b4"  ))

  }

}
