#ifndef KLSAMPLESCOST_H_
#define KLSAMPLESCOST_H_

#include "MUQ/Optimization/CostFunction.h"

namespace muq {
  namespace Approximation {
    class KLSamplesCost : public muq::Optimization::CostFunction {
    public:

      KLSamplesCost(Eigen::MatrixXd const& vand, Eigen::MatrixXd const& deriv);

      virtual ~KLSamplesCost() = default;

      virtual void SetPoint(Eigen::VectorXd const& evalPt) override;
      virtual double Cost() override;

      virtual Eigen::VectorXd Gradient() override;
      virtual Eigen::MatrixXd Hessian() override;
      virtual Eigen::VectorXd ApplyHessian(Eigen::VectorXd const& vec) override;

    private:

      Eigen::MatrixXd const& vand;
      Eigen::MatrixXd const& deriv;

      const double numSamps;

      Eigen::VectorXd vandApp; // holds V*x
      Eigen::VectorXd derivApp; // holds G*x

      Eigen::MatrixXd vv; // holds vand.transpose()*vand

    };

    class KLSamplesConstraint : public muq::Modeling::ModPiece {
      public:
      KLSamplesConstraint(Eigen::MatrixXd const& deriv);

      virtual ~KLSamplesConstraint() = default;

     private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

      virtual void JacobianImpl(unsigned int const outputDimWrt, unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

      const Eigen::MatrixXd& deriv;
    };
  } // namespace Approximation
} // namespace muq

#endif
