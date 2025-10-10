#include "cardillo.hpp"

namespace cardillo {

class PointMass {
    private:
        real_t mass_;
        Matrix33r M_;
        Matrix33r M_inv_;

    public:
        PointMass(real_t mass) :  mass_{mass} {
            M_.diagonal() << mass, mass, mass;
            M_inv_.diagonal() << 1 / mass, 1 / mass, 1 / mass;
        };
};

}  // namespace cardillo