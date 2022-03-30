#include "MUQ/Modeling/HTTPModel/HTTPModPiece.h"
#include "MUQ/Modeling/HTTPModel/HTTPComm.h"

namespace muq {
  namespace Modeling {

    /**
      @class HTTPModPieceWrapper
      @brief Wrap a ModPiece in an UM-Bridge Model
      @details This is needed in order to easily serve a MUQ ModPiece via UM-Bridge.
      */
    class HTTPModPieceWrapper : public ShallowModPiece {
    public:

      HTTPModPieceWrapper(std::shared_ptr<muq::Modeling::ModPiece> modPiece)
      : ShallowModPiece(modPiece->inputSizes, modPiece->outputSizes), modPiece(modPiece)
      {}

      void Evaluate(std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& inputs, json config) override {
        outputs = modPiece->Evaluate(inputs);
      }

      void Gradient(unsigned int outWrt,
                    unsigned int inWrt,
                    std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& inputs,
                    Eigen::VectorXd const& sens,
                    json config = json()) override {
        gradient = modPiece->Gradient(outWrt, inWrt, inputs, sens);

      }

      void ApplyJacobian(unsigned int outWrt,
                    unsigned int inWrt,
                    std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& inputs,
                    Eigen::VectorXd const& vec,
                    json config = json()) override {
        jacobianAction = modPiece->ApplyJacobian(outWrt, inWrt, inputs, vec);
      }

      void ApplyHessian(unsigned int outWrt,
                        unsigned int inWrt1,
                        unsigned int inWrt2,
                        muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs,
                        Eigen::VectorXd const& sens,
                        Eigen::VectorXd const& vec,
                        json config = json()) override {
        hessAction = modPiece->ApplyHessian(outWrt, inWrt1, inWrt2, inputs, sens, vec);
      }

      bool SupportsEvaluate() override {return true;}
      bool SupportsGradient() override {return true;} // Expose derivative information as well,
      bool SupportsApplyJacobian() override {return true;} // since we fall back to finite differences automatically
      bool SupportsApplyHessian() override {return true;}

    private:
      std::shared_ptr<muq::Modeling::ModPiece> modPiece;
    };

    /**
     * @brief Serve a ModPiece via network using UM-Bridge
     *
     * @param modPiece The modPiece to serve via UM-Bridge
     * @param host Bind address, may be 0.0.0.0
     * @param port Port at which to serve the modPiece
     */
    void serveModPiece(std::shared_ptr<ModPiece> modPiece, std::string host, int port) {
      HTTPModPieceWrapper wrapper(modPiece);
      ::serveModPiece(wrapper, host, port);
    }

  }
}
