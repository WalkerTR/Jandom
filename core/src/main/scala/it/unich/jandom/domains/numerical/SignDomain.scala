package it.unich.jandom.domains.numerical

import it.unich.jandom.domains.CachedTopBottom
import it.unich.jandom.domains.WideningDescription
// import it.unich.jandom.utils.numberext.RationalExt
// import spire.math.Rational

/**
 * This is the domain of signs.
 *
 * @constructor Builds a box domain using doubles as bounds
 * @param overReals is true if the domain is correct w.r.t. real arithmetic, otherwise it is correct w.r.t.
 * double arithmetics.
 */
class SignDomain extends NumericalDomain {

  /**
    * This is the ADT representing abstract elements of the domain, i.e. signs
    */
  sealed trait Sign extends PartiallyOrdered[Sign]
  case class TopSign() extends Sign {
    def tryCompareTo[B >: SignDomain.this.Sign](that : B)
        (implicit evidence$1: B => scala.math.PartiallyOrdered[B]): Option[Int] = {
      that match {
        case other: Sign => other match {
          case TopSign() => Option(0)
          case _ => Option(1)
        }
        case _ => Option.empty
      }
    }
  }
  case class PlusSign() extends Sign {
    def tryCompareTo[B >: SignDomain.this.Sign](that : B)
        (implicit evidence$1: B => scala.math.PartiallyOrdered[B]): Option[Int] = {
      that match {
        case other: Sign => other match {
          case PlusSign() => Option(0)
          case TopSign() => Option(-1)
          case BotSign() => Option(1)
          case _ => Option.empty
        }
        case _ => Option.empty
      }
    }
  }
  case class MinusSign() extends Sign {
    def tryCompareTo[B >: SignDomain.this.Sign](that : B)
        (implicit evidence$1: B => scala.math.PartiallyOrdered[B]): Option[Int] = {
      that match {
        case other: Sign => other match {
          case MinusSign() => Option(0)
          case TopSign() => Option(-1)
          case BotSign() => Option(1)
          case _ => Option.empty
        }
        case _ => Option.empty
      }
    }
  }
  case class ZeroSign() extends Sign {
    def tryCompareTo[B >: SignDomain.this.Sign](that : B)
        (implicit evidence$1: B => scala.math.PartiallyOrdered[B]): Option[Int] = {
      that match {
        case other: Sign => other match {
          case ZeroSign() => Option(0)
          case TopSign() => Option(-1)
          case BotSign() => Option(1)
          case _ => Option.empty
        }
        case _ => Option.empty
      }
    }
  }
  case class BotSign() extends Sign {
    def tryCompareTo[B >: SignDomain.this.Sign](that : B)
        (implicit evidence$1: B => scala.math.PartiallyOrdered[B]): Option[Int] = {
      that match {
        case other: Sign => other match {
          case BotSign() => Option(0)
          case _ => Option(-1)
        }
        case _ => Option.empty
      }
    }
  }

  private def signUnion (x : Sign, y : Sign) : Sign =
    x.tryCompareTo(y) match {
      case Some(0) => x
      case Some(-1) => y
      case Some(1) => x
      case None => TopSign()
      case _ => throw new IllegalArgumentException()
    }

  private def signIntersection (x : Sign, y : Sign) : Sign =
    x.tryCompareTo(y) match {
      case Some(0) => x
      case Some(-1) => x
      case Some(1) => y
      case None => BotSign()
      case _ => throw new IllegalArgumentException()
    }

  private def signToInterval (x : Sign) : (Double, Double) =
    x match {
      case TopSign() => (Double.NegativeInfinity, Double.PositiveInfinity)
      case PlusSign() => (Double.MinPositiveValue, Double.PositiveInfinity)
      case MinusSign() => (Double.NegativeInfinity, -Double.MinPositiveValue)
      case ZeroSign() => (0,0)
      case BotSign() => (Double.PositiveInfinity, Double.NegativeInfinity)
    }

  private def intervalToSign (x: Double, y: Double): Sign =
    if (x == 0 && y == 0) ZeroSign() else
      if (x < 0 && 0 < y) TopSign() else
        if (x < 0 && y < 0) MinusSign() else
          if (x > 0 && y > 0) PlusSign() else
            BotSign()

  /**
   * This is the class representing a single abstract element, i.e. a sign.
   *
   * @constructor Creates a sign.
   * @param sign the sign.
   * @param high the upper bounds of the box.
   * @param isEmpty is true when the box is empty. It is needed for the case of 0-dimensional boxes.
   * @note `low`, `high` and `isEmpty` should be normalized according to the method `normalized`
   * @throws IllegalArgumentException if parameters are not correct.
   */
  final class Property(val sign: Seq[Sign])
      extends NumericalProperty[Property] {

    type Domain = SignDomain
    def domain = SignDomain.this

    def indexWithinDimension(n : Int) : Boolean = (0 <= n) && (n < dimension)

    private def asInterval = {
      val (x1, x2) = this.sign.map(signToInterval).toArray.unzip
      BoxDoubleDomain(false).apply(x1, x2)
    }

    def fromInterval (interval : BoxDoubleDomain#Property) : Property =
      new Property(
        (interval.low, interval.high).zipped.map(intervalToSign))

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def union(that: Property): Property =
      new Property((this.sign, that.sign).zipped map signUnion)

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def intersection(that: Property): Property =
      new Property((this.sign, that.sign).zipped map signIntersection)

    /**
     * This is the standard widening on signs, i.e. one that does nothing, since
     * the domain is finite.
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def widening(that: Property) = that union this

    /**
     * This is the standard narrowing on signs, i.e. one that does nothing,
     * since the domain is finite.
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def narrowing(that: Property) = that

    def minimize(lf: LinearForm) = ??? // RationalExt(linearEvaluation(lf)._1)

    def maximize(lf: LinearForm) = ??? // RationalExt(linearEvaluation(lf)._2)

    def frequency(lf: LinearForm) = ???

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def nonDeterministicAssignment(n: Int): Property =
      fromInterval(asInterval.nonDeterministicAssignment(n))

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @todo @inheritdoc
     * @throws $ILLEGAL
     */
    def linearAssignment(n: Int, lf: LinearForm): Property =
      fromInterval(asInterval.linearAssignment(n, lf))

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @todo @inheritdoc
     * @throws $ILLEGAL
     */
    def linearInequality(lf: LinearForm): Property =
      fromInterval(asInterval.linearInequality(lf))

    def constraints = asInterval.constraints

    def isPolyhedral = true

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def linearDisequality(lf: LinearForm): Property =
      fromInterval(asInterval.linearDisequality(lf))

    /**
     * @inheritdoc
     * This is a complete operator for boxes.
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def addVariable: Property = SignDomain.this(this.sign :+ TopSign())

    def removeAt[A](i : Int, arr : Seq[A]) : Seq[A] =
      arr.zipWithIndex.filter(_._2 != i).map(_._1)

    /**
     * @inheritdoc
     * This is a complete operator for boxes.
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def delVariable(n: Int): Property = {
      require(n < dimension && n >= 0)
      new Property(removeAt(n, this.sign))
    }

    /**
     * @inheritdoc
     * This is a complete operator for boxes.
     * @note @inheritdoc
     * @throws IllegalArgumentException if parameters are not correct (but we do not check injectivity of `rho`)
     */
    def mapVariables(rho: Seq[Int]) = {
      require(rho.length == dimension)
      val newdim = rho.count(_ >= 0)
      require(rho forall { i => i >= -1 && i < newdim })
      // we do not check injectivity
      val newSign = new Array[Sign](newdim)
      for ((newi, i) <- rho.zipWithIndex; if newi >= 0) {
        newSign(newi) = sign(i)
      }
      new Property(newSign)
    }

    /**
     * @inheritdoc
     * @throws $ILLEGAL
     */
    def mkString(vars: Seq[String]): String = {
      require(vars.length >= dimension)
      if (isEmpty)
        "empty"
      else {
        val bounds = for (i <- 0 until dimension) yield {
          sign(i) match {
            case TopSign() => s"${vars(i)} don't know"
            case PlusSign() => s"${vars(i)} positive"
            case MinusSign() => s"${vars(i)} negative"
            case ZeroSign() => s"${vars(i)} zero"
            case BotSign() => s"${vars(i)} unreachable"
          }
        }
        bounds.mkString("[ ", " , ", " ]")
      }
    }

    val dimension: Int = sign.length


    def isEmpty = isBottom

    def isBottom = this.sign.forall(_== BotSign())
    def isTop = this.sign.forall(_== TopSign())

    def bottom = SignDomain.this.bottom(sign.length)
    def top = SignDomain.this.top(sign.length)

    def tryCompareTo[B >: Property](other: B)(
      implicit arg0: (B) => PartiallyOrdered[B]): Option[Int] =

      other match {
        case other: Property =>
          require(dimension == other.dimension)
          (isEmpty, other.isEmpty) match {
            case (true, true) => Option(0)
            case (false, true) => Option(1)
            case (true, false) => Option(-1)
            case (false, false) =>
              val pairs = (this.sign, other.sign).zipped
              if (pairs.forall(_ == _))
                Option(0)
              else if (pairs.forall(_ <= _))
                Option(-1)
              else if (pairs.forall(_ >= _))
                Option(1)
              else
                Option.empty
          }
        case _ => Option.empty
    }

    override def hashCode: Int = 41 * (41 + sign.hashCode) + sign.hashCode
  }

  val widenings = Seq(WideningDescription.default[Property])

  /**
   * Returns a sign.
   * @return the abstract sign with the specified signs.
   * @throws $ILLEGAL
   */
  def apply(sign: Seq[Sign]): Property = new Property(sign)

  private def signOf(d: Double): Sign =
    if (d == 0) ZeroSign() else
      if (d < 0) MinusSign() else
        PlusSign()

  // /**
  //  * Returns a box consisting of the single point `poDouble`.
  //  */
  // def apply(poDouble: Seq[Double]): Property = apply(poDouble.map(signOf))

  /**
   * @inheritdoc
   * @note @inheritdoc
   * @throws $ILLEGAL
   */
  def top(n: Int): Property = new Property(Array.fill(n)(TopSign()))

  /**
   * @inheritdoc
   * @note @inheritdoc
   * @throws $ILLEGAL
   */
  def bottom(n: Int): Property = new Property(Array.fill(n)(BotSign()))
}

object SignDomain {
  /**
   * Returns an abstract domain for boxes which is correct w.r.t. real arithmetic or
   * double arithmetic, according to the parameter `overReals`.
   */
  def apply() = this.overDoubles

  /**
   * The domain of boxes correct w.r.t. doubles and with cached top and bottom.
   */
  private val overDoubles = new SignDomain() with CachedTopBottom
}
