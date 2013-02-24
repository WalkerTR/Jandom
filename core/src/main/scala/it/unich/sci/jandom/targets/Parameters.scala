/**
 * Copyright 2013 Gianluca Amato
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

package it.unich.sci.jandom
package targets
  
import domains.{NumericalProperty,NumericalDomain}
import widenings.{Widening, DefaultWidening}
import narrowings.{Narrowing, DefaultNarrowing}
import ppfactories.PPFactory
import it.unich.sci.jandom.parameters.WideningScope
import it.unich.sci.jandom.parameters.NarrowingStrategy

/**
 * This class is used to keep parameters for analyzers.
 * @tparam Tgt the target related to this parameter
 * @param domain the numerical domain for the analysis
 * @param tgt the target for the analysis
 * @author Gianluca Amato <amato@sci.unich.it>
 *
 */
class Parameters[Tgt <: Target] (val domain: NumericalDomain, val tgt: Tgt) {
  
  type Property = domain.Property
  
  /**
  * The widening factory used in the analysis. Defaults to the factory for the standard domain widening.
  */
  var wideningFactory: PPFactory[Tgt, Widening] = DefaultWidening
  
  /**
   * The narrowing factory used in the analysis. Defaults to the standard domain narrowing. 
   */
  var narrowingFactory: PPFactory[Tgt,Narrowing] = DefaultNarrowing
  
  /**
   * This parameter determines whether results are saved for each program point or only for widening points.
   */
  var allPPResult = true
  
  /**
   * This parameter determines whether standard or local widening is used. At the moment, this is only supported
   * by the SLSL target.
   */
  var wideningScope = WideningScope.Output
  
  /**
   * This parameter determine the interlacing strategy between narrowing and widening
   */
  var narrowingStrategy = NarrowingStrategy.Restart
  
  /**
   * This is used for putting results in tags
   */
  var tag = scala.collection.mutable.Map[Any, Property]()
  
  /** 
   * This is a variable globally used by the analyzer for keeping track of nested level
   */
  var nestingLevel = 0
  
  /**
   * This is a java writer where the analyzer write debug informations
   */
  var debugWriter = new java.io.Writer {
    override def write(cbuf: Array[Char], off: Int, len: Int) {}
    override def flush() { }
    override def close() { }
    override def toString = ""
  }
  
  def log(msg: String) {
    debugWriter.write(" "*nestingLevel*3) 
    debugWriter.write(msg)
  }
}
