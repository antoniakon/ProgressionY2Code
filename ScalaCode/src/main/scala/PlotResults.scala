

/**
  * Created by Antonia Kontaratou.
  */
class PlotResults {
  // Function to plot the results of the Gibbs sampler (based on Scala for data science book, p. 63)
  // Plot the results
  def plotResults(resultsMatrix:DenseMatrix[Double], labels:List[String]): Unit = {
    val time = linspace(1, resultsMatrix.rows, resultsMatrix.rows) // linspace creates a denseVector of Doubles. Otherwise time = convert(DenseVector.range(1,st.rows+1,1), Double), bcs if it is not converted to Double the plot below is not created. Error for yv implisit.
    val fig = Figure("Gibbs sampler results")
    val ncols = resultsMatrix.cols
    require(ncols == labels.size, "Number of columns in feature matrix must match length of labels.")
    fig.clear


    (0 until ncols).foreach { irow =>
      val pTrace = fig.subplot(3, ncols, irow)
      tracePlot(pTrace)(time, resultsMatrix(::, irow), "time", labels(irow))
      val pHist= fig.subplot(3, ncols,irow+ncols)
      plotHistogram(pHist)(resultsMatrix(::,irow), labels(irow))
      val pAutocor= fig.subplot(3, ncols,irow+2*ncols)
      plotAutocorrelation(pAutocor)(resultsMatrix(::,irow), labels(irow))
    }

    // Function for the traceplots
    def tracePlot(plt: Plot)(xdata: DenseVector[Double], ydata:DenseVector[Double], xlabel:String, ylabel: String): Unit ={
      plt += plot(xdata, ydata, '-')
      plt.xlabel = xlabel
      plt.ylabel = ylabel
      plt.title = "Trace of " + ylabel
    }
    // Function for the histograms
    def plotHistogram(plt:Plot)(data: DenseVector[Double], label:String): Unit ={
      plt += hist(data, bins=20)
      plt.xlabel = label
      plt.title = "Histogram of "+label
    }

    // Function for the autocorrelation plots
    def plotAutocorrelation(plt:Plot)(data: DenseVector[Double], ylabel:String): Unit ={
      val maxLag=25
      val lag = linspace(0, maxLag, maxLag)
      val dataSize = data.length

      def estimateACF(data: DenseVector[Double], maxLag: Int): DenseVector[Double] ={
        val acfRes= new DenseVector[Double](maxLag)
        var x = new DenseVector[Double](data.length)
        x= data-mean(data)
        //acfRes(0)=sum(x(0 to dataSize-1):*x(0 to dataSize-1)) / sum(pow(x,2))
        (0 until maxLag).foreach{K=>
          val M = dataSize-K-1
          acfRes(K)= sum(x(0 to M):*x((K) to (M+K))) / sum(pow(x,2))
        }
        acfRes
      }

      val acfData= estimateACF(data, maxLag)
      plt += plot(lag, acfData, '.')
      plt.xlabel = "lag"
      plt.ylabel = ylabel
      plt.title = "Autocorrelation plot of "+ylabel
      //The correlation coefficient for the scatterplot summarizes the strength of the linear relationship between present and past values. The correlation coefficient is annotated at the top of the plot, along with the correlation that would be considered approximately significant at the 95% significance level against the null hypothesis that true correlation is zero.
      // if a series is completely random, then, for large sample size sample size N, the lagged-correlation coefficient is approximately normally distributed with mean 0 and variance 1/N. The probability is thus roughly ninety-five percent that the correlation falls within two standard deviations, or 2.0/sqrt(N
      val upconf= 1.96/ sqrt(dataSize) // For a 95%-confidence interval, the critical value is  and the confidence interval is sqrt(2)~1.96
      val lconf= -1.96/ sqrt(dataSize)
      plt.plot.addRangeMarker(new ValueMarker(upconf))
      plt.plot.addRangeMarker(new ValueMarker(lconf))
      plt.plot.getRangeAxis.setUpperBound(1.0)
      plt.plot.getRangeAxis.setLowerBound(-1.0)

    }
  }
}
