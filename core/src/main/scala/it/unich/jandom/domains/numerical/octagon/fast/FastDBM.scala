package it.unich.jandom.domains.numerical.octagon.fast

import breeze.numerics.pow
import it.unich.jandom.domains.numerical.octagon._

import scala.language.higherKinds
import scalaz.{Applicative, Apply, Monoid, Traverse}
import scalaz.std.option._
import scalaz.std.list._

// Sparsity D = 1 - (nni / (2n^2 + 2n))
// Switching between DBMs: the sparsity can increase, for instance during
// widening. Recovering sparsity information and  independent components has a
// quadratic worst case complexity, so we only perform it by piggybacking on
// the closure operator. We also use closure computations as switching points.


// ADT of "closable" DBMs in their fast implementation from Vechev et al.
// They are "closable" in the sense that they augment the ADT of fast DBMs with
// the type-level capability of being indexed by their strong closure state.


class FastDBM[A](private val m: HalfMatrix[A],
                 private val components: Set[Set[VarIndex]],
                 val dimension: Int)
                (implicit ifield: InfField[A]) {

  import VarIndexOps._

  val vars: Set[VarIndex] = components.foldRight[Set[VarIndex]](Set())(_ ++ _)

  def update(f: (Int, Int) => A): FastDBM[A] = {
    val varIndeces = for (i <- 0 until dimension) yield VarIndex(i)
    new FastDBM(m.update(f), Set(varIndeces.toSet), dimension)
  }

  def get(i: Int, j: Int): A =
    if (vars.contains(VarIndex(i/2)) && vars.contains(VarIndex(j/2))) {
      m(i, j)
    } else if (i == j) {
      ifield.zero
    } else {
      ifield.infinity
    }

  def strongClosure: FastDBM[A] = ???
  def incrementalClosure(v: VarIndex): FastDBM[A] = ???

  def forget(vi: VarIndex): FastDBM[A] = {
    val newComponents = components.map(_.filter(_ != vi))
    new FastDBM(m, newComponents, dimension)
  }

  def flipVar(vi: VarIndex): FastDBM[A] =
    if (vars.contains(vi)) {

      var updateVars: Set[(VarIndex, VarIndex)] =
        (vars.map((_, vi)) union (vars.map((vi, _))))
        .filter({ case (v1, v2) => v1 <= v2 })

      val indeces = varsToIndeces(updateVars)

      val f: (Int, Int) => A = (i, j) => {
        if (i == varPlus(vi) || i == varMinus(vi)) {
          if (j == varPlus(vi) || j == varMinus(vi))
            m.get(signed(i), signed(j))
          else
            m.get(signed(i), j)
        } else {
          if (j == varPlus(vi) || j == varMinus(vi))
            m.get(i, signed(j))
          else
            m.get(i, j)
        }
      }

      // apply f on m_iv and m_vj elements only, where
      // i and j (and v) are used
      val newM = indeces.foldLeft(m)((newM, idx) => {
        val (i, j) = idx
        newM.update(i, j, f(i, j))
      })

      new FastDBM(newM, components, dimension)
    } else {
      this
    }

  def addScalarOnVar(vi: VarIndex, const: A): FastDBM[A] = 
    if (vars.contains(vi)) {

      var updateVars: Set[(VarIndex, VarIndex)] =
        (vars.map((_, vi)) union (vars.map((vi, _))))
        .filter({ case (v1, v2) => v1 <= v2 })

      val indeces = varsToIndeces(updateVars)

      val f: (Int, Int) => A = (i, j) => {
        val g1 = (i == varPlus(vi) && j != varPlus(vi) && j != varMinus(vi)) ||
                 (j == varMinus(vi) && i != varPlus(vi) && i != varMinus(vi))
        val g2 = (i != varPlus(vi) && i != varMinus(vi) && j == varPlus(vi)) ||
                 (j != varPlus(vi) && j != varMinus(vi) && i == varMinus(vi))
        val g3 = i == varPlus(vi) && j == varMinus(vi)
        val g4 = i == varMinus(vi) && j == varPlus(vi)
        if (g1) ifield.-(m.get(i, j), const) else
        if (g2) ifield.+(m.get(i, j), const) else
        if (g3) ifield.-(m.get(i, j), ifield.double(const)) else
        if (g4) ifield.+(m.get(i, j), ifield.double(const)) else
          m.get(i, j)
      }

      // apply f on m_iv and m_vj elements only, where
      // i and j (and v) are used
      val newM = indeces.foldLeft(m)((newM, idx) => {
        val (i, j) = idx
        newM.update(i, j, f(i, j))
      })

      new FastDBM(newM, components, dimension)
    } else {
      this
    }


  def dbmUnion(other: FastDBM[A]): FastDBM[A] = {
    val newVars = vars intersect other.vars
    val newComponents = components.map(_.filter(newVars.contains(_)))

    val updateVars = for(i <- newVars; j <- newVars) yield (i, j)
    val indeces = varsToIndeces(updateVars)

    val f: (Int, Int) => A = (i, j) =>
      ifield.max(m.get(i, j), other.m.get(i, j))

    // apply f on m_iv and m_vj elements only, where
    // i and j (and v) are used
    val newM = indeces.foldLeft(m)((newM, idx) => {
      val (i, j) = idx
      newM.update(i, j, f(i, j))
    })

    new FastDBM(newM, newComponents, dimension)
  }

  def dbmIntersection(other: FastDBM[A]): FastDBM[A] = {
    val newVars = vars union other.vars
    val newComponents = Set(newVars)
    
    val updateVars = for(i <- newVars; j <- newVars) yield (i, j)
    val indeces = varsToIndeces(updateVars)

    val f: (Int, Int) => A = (i, j) =>
      ifield.min(m.get(i, j), other.m.get(i, j))

    // apply f on m_iv and m_vj elements only, where
    // i and j (and v) are used
    val newM = indeces.foldLeft(m)((newM, idx) => {
      val (i, j) = idx
      newM.update(i, j, f(i, j))
    })

    new FastDBM(newM, newComponents, dimension)
  }

  def widening(other: FastDBM[A]): FastDBM[A] = ???
  def narrowing(other: FastDBM[A]): FastDBM[A] = ???

  def isTop: Boolean = {
    val checkVars = for(i <- vars; j <- vars; if i <= j) yield (i, j)
    val indeces = varsToIndeces(checkVars)
    indeces.forall({ case (i, j) =>
      if (i == j) {
        get(i, j) == ifield.zero
      } else {
        get(i, j) == ifield.infinity
      }
    })
  }
  def isBottom: Boolean = ???


  // convert couples of variable indeces to couples of element indeces
  // e.g. {(0,1)} to {(0,2), (0,3), (1,2), (1,3)}
  private def varsToIndeces(variables: Set[(VarIndex, VarIndex)]) =
    variables.flatMap({ case (v1, v2) =>
      Set(
        (varPlus(v1), varPlus(v2)),
        (varPlus(v1), varMinus(v2)),
        (varMinus(v1), varPlus(v2)),
        (varMinus(v1), varMinus(v2))
      )
    })

}