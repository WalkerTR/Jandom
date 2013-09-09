/**
 * Copyright 2013 Gianluca Amato <gamato@unich.it>
 *
 * This file is part of JANDOM: JVM-based Analyzer for Numerical DOMains
 * JANDOM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * JANDOM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty ofa
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with JANDOM.  If not, see <http://www.gnu.org/licenses/>.
 */

package it.unich.sci.jandom.domains.numerical

import scala.collection.mutable.WeakHashMap

/**
 * The most abstract numerical domain. Not that it does not extends from [[it.unich.sci.jandom.domains.TopLike]]
 * since numerical domains do not have a well defined top and bottom domains, but are fibered over dimensions.
 * @author Gianluca Amato <gamato@unich.it>
 */
object NumericalTopDomain extends NumericalDomain {

  /**
   * This is the class for the top property. There is a single top property for each dimension.
   */
  case class Property private[NumericalTopDomain] (val dimension: Int) extends NumericalProperty[Property] {
    def tryCompareTo[B >: Property](other: B)(implicit arg0: (B) => PartiallyOrdered[B]): Option[Int] = other match {
      case other: Property => Some(0)
      case _ => None
    }
    def widening(that: Property) = { require(that.dimension == dimension); this }
    def union(that: Property) = { require(that.dimension == dimension); this }
    def narrowing(that: Property) = { require(that.dimension == dimension); this }
    def intersection(that: Property) = { require(that.dimension == dimension); this }
    def nonDeterministicAssignment(n: Int) = { require(n < dimension); this }
    def linearAssignment(n: Int, lf: LinearForm[Double]) = {
      require(n < dimension && lf.dimension <= dimension)
      this
    }
    def linearInequality(lf: LinearForm[Double]) = { require(lf.dimension <= dimension); this }
    def linearDisequality(lf: LinearForm[Double]) = { require(lf.dimension <= dimension); this }

    def minimize(lf: LinearForm[Double]) = {
      require(lf.dimension <= dimension)
      if (lf.homcoeffs.exists(_ != 0))
        Double.NegativeInfinity
      else
        lf.known
    }

    def maximize(lf: LinearForm[Double]) = {
      require(lf.dimension <= dimension)
      if (lf.homcoeffs.exists(_ != 0))
        Double.PositiveInfinity
      else
        lf.known
    }

    def frequency(lf: LinearForm[Double]) = {
      require(lf.dimension <= dimension)
      if (lf.homcoeffs.exists(_ != 0))
        None
      else
        Some(lf.known)
    }

    def addVariable = new Property(dimension + 1)
    def delVariable(n: Int) = { require(n < dimension); new Property(dimension - 1) }
    def mapVariables(rho: Seq[Int]) = new Property(rho count { _ != -1 })
    def isEmpty = false
    def isTop = true
    def isBottom = true
    def bottom = this
    def top = this
    def mkString(vars: Seq[String]): String = "top"
  }

  def top(n: Int) = Property(n)
  def bottom(n: Int) = Property(n)
}