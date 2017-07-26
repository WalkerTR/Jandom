package it.unich.jandom.domains.numerical.octagon.fast

/**
  * Created by fsestini on 7/11/17.
  *
  * Implementation of matrices that are represented by storing only the lower
  * triangular half of it, as explained in Singh et al.
  */
class HalfMatrix[A](private val vec: Vector[A], val dimension: Int) {

  require(dimension > 0)

  private def getIndex(i: Int, j: Int): Int = j + ((i+1)*(i+1))/2

  private def elementIndex(i: Int, j: Int): Int =
    if (i < j)
      getIndex(j^1, i^1)
    else
      getIndex(i, j)

  def get(i: Int, j: Int): A = vec(elementIndex(i,j))

  def update(i: Int, j: Int, x: A): HalfMatrix[A] = {
    require(0 <= i && i < dimension && 0 <= j && j < dimension)
    new HalfMatrix(vec.updated(elementIndex(i, j), x), dimension)
  }

  def update(updater: (Int, Int) => A): HalfMatrix[A] = {
    val values = for (i <- 0 until dimension;
                      j <- 0 to (if (i % 2 == 0) i+1 else i))
                  yield updater(i, j)
    new HalfMatrix(values.toVector, dimension)
  }

  def apply(i: Int, j: Int): A = {
    require(0 <= i && i < dimension && 0 <= j && j < dimension)
    vec(elementIndex(i, j))
  }

  def combine[B, C](that: HalfMatrix[B], f: (A, B) => C): HalfMatrix[C] = {
    require(this.dimension == that.dimension)
    val mat = (this.vec zip that.vec) map (f.tupled)
    new HalfMatrix(mat, dimension)
  }

}
