// Range of +Grid inter-plane links over one orbit from the toolkit's library propagator (Kepler + J2), unmodified.
// Prints JSON lines: shell a b [u_deg...] [range_m...]. No routing, no packets; used only to validate Lemma 1 (m10).
#include "ns3/core-module.h"
#include "ns3/sgp4-mobility-model.h"
#include "ns3/walker-constellation.h"
#include <cmath>
#include <iostream>
using namespace ns3;
using namespace ns3::ntncon;
int main()
{
    struct S { const char* n; uint32_t P, Q, F; double i, h; bool star; };
    for (const S& sh : {S{"starlink", 72, 22, 39, 53.0, 550.0, false}, S{"kuiper", 34, 34, 1, 51.9, 630.0, false},
                        S{"oneweb", 12, 49, 1, 87.9, 1200.0, true}})
    {
        WalkerConfig wc;
        wc.num_planes = sh.P; wc.total_sats = sh.P * sh.Q; wc.phasing_f = sh.F; wc.inclination_deg = sh.i;
        wc.altitude_km = sh.h; wc.epoch_unix_s = 1792800000.0;
        auto els = sh.star ? WalkerConstellation::BuildStar(wc) : WalkerConstellation::BuildDelta(wc);
        std::vector<Ptr<Sgp4MobilityModel>> m;
        for (auto& e : els)
        {
            auto x = CreateObject<Sgp4MobilityModel>();
            x->SetElements(e); x->SetSimulationEpochUnix(wc.epoch_unix_s);
            m.push_back(x);
        }
        const double r = 6378137.0 + sh.h * 1e3;
        const double T = 2 * M_PI * std::sqrt(r * r * r / 3.986004418e14);
        for (uint32_t s : {0u, sh.Q / 3, 2 * sh.Q / 3})
        {
            const uint32_t a = s, b = sh.Q + s;
            std::cout << "{\"shell\":\"" << sh.n << "\",\"a\":" << a << ",\"b\":" << b << ",\"T\":" << T << ",\"range\":[";
            for (int k = 0; k <= 200; ++k)
            {
                Vector pa, pb, v;
                m[a]->GetEcefStateAt(Seconds(T * k / 200.0), pa, v);
                m[b]->GetEcefStateAt(Seconds(T * k / 200.0), pb, v);
                const double d = std::sqrt(std::pow(pa.x - pb.x, 2) + std::pow(pa.y - pb.y, 2) + std::pow(pa.z - pb.z, 2));
                std::cout << (k ? "," : "") << std::fixed << d;
            }
            std::cout << "]}\n";
        }
    }
    return 0;
}
