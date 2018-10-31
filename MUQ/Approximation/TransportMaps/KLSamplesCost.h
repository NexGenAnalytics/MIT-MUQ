#ifndef KLSAMPLESCOST_H_
#define KLSAMPLESCOST_H_

#include "MUQ/Optimization/CostFunction.h"

namespace muq {
  namespace Approximation {
    class KLSamplesCost : public muq::Optimization::CostFunction {
    public:

      KLSamplesCost(Eigen::MatrixXd const& vand, Eigen::MatrixXd const& deriv);

      virtual ~KLSamplesCost() = default;

    private:

      virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

      void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override;

      const Eigen::MatrixXd& vand;

      const Eigen::MatrixXd& deriv;
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
