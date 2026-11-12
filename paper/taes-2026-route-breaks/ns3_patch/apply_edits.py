#!/usr/bin/env python3
"""Produce the patched scheduler and router from orig/ (the example is edited by hand; see README)."""
import re, shutil, os
H = os.path.dirname(os.path.abspath(__file__))
def R(f): return open(os.path.join(H, "orig", f)).read()
def W(f, s): open(os.path.join(H, f), "w").write(s)
def rep(s, a, b, n=1):
    assert s.count(a) >= 1 and (n == 0 or s.count(a) == n), a
    return s.replace(a, b)

# ---------- scheduler
h = R("contact-graph-scheduler.h")
h = rep(h, "    /// Start the periodic visibility sampler.", """    /// P7: the ground links that can occur in the run. With a list set, only these are sampled, and a
    /// satellite well below the mask is skipped until it could first reach it at the rate bound below;
    /// transitions still land on the same sampling grid.
    void SetGslCandidates(const std::vector<std::pair<uint32_t, uint32_t>>& satGs);
    /// P7: refresh the candidate ISLs only every k-th tick (k = 1 is the original behavior).
    void SetIslRefreshEvery(uint32_t k) { m_islEvery = k; }

    /// Start the periodic visibility sampler.""")
h = rep(h, "    uint64_t m_gslUp{0};", """    struct GslCand
    {
        uint32_t sat;
        uint32_t gs;
        Ptr<Sgp4MobilityModel> mob;
        double latDeg;
        double lonDeg;
        Vector gsEcef;
        bool up{false};
        double nextCheck{0.0};
    };
    std::vector<GslCand> m_gslCand;
    uint32_t m_islEvery{1};
    uint64_t m_tickN{0};

    uint64_t m_gslUp{0};""")
W("contact-graph-scheduler.h", h)
c = R("contact-graph-scheduler.cc")
c = rep(c, 'NS_LOG_COMPONENT_DEFINE("ContactGraphScheduler");', '''#ifdef P7_STANDALONE
NS_LOG_COMPONENT_DEFINE("ContactGraphSchedulerP7");
#else
NS_LOG_COMPONENT_DEFINE("ContactGraphScheduler");
#endif''')
c = rep(c, 'TypeId("ns3::ntncon::ContactGraphScheduler")', '''TypeId(
#ifdef P7_STANDALONE
                                   "ns3::ntncon::ContactGraphSchedulerP7"
#else
                                   "ns3::ntncon::ContactGraphScheduler"
#endif
                                   )''')
c = rep(c, """void
ContactGraphScheduler::Start()
{""", """void
ContactGraphScheduler::SetGslCandidates(const std::vector<std::pair<uint32_t, uint32_t>>& satGs)
{
    m_gslCand.clear();
    for (const auto& [s, gsId] : satGs)
    {
        const auto& g = m_gss.at(gsId);
        const double lat = g.lat_deg * M_PI / 180.0;
        const double lon = g.lon_deg * M_PI / 180.0;
        const double Ngs = kEarthRadiusM / std::sqrt(1.0 - kWgs84E2 * std::sin(lat) * std::sin(lat));
        GslCand gc;
        gc.sat = s;
        gc.gs = gsId;
        gc.mob = m_sats.at(s);
        gc.latDeg = g.lat_deg;
        gc.lonDeg = g.lon_deg;
        gc.gsEcef = Vector(Ngs * std::cos(lat) * std::cos(lon),
                           Ngs * std::cos(lat) * std::sin(lon),
                           Ngs * (1.0 - kWgs84E2) * std::sin(lat));
        m_gslCand.push_back(gc);
    }
}

void
ContactGraphScheduler::Start()
{""")
# GSL loop: replace the whole block when a candidate list is set
c = rep(c, """    // GSL visibility.
    for (auto& gsIt : m_gss)
    {""", """    // P7: sparse GSL sampling. Elevation changes no faster than kMaxElevRate deg/s at any LEO altitude
    // used here (0.9 deg/s overhead at 460 km), so a satellite e degrees below the mask cannot cross it
    // for (e / rate) seconds. Decisions stay on the sampling grid and match the full scan.
    constexpr double kMaxElevRate = 2.0;
    for (auto& c : m_gslCand)
    {
        if (!c.up && t < c.nextCheck)
        {
            continue;
        }
        const double elev = c.mob->GetElevationDeg(c.latDeg, c.lonDeg);
        const bool visible = c.up ? (elev >= m_minElevDeg - m_gateHysteresisDeg) : (elev >= m_minElevDeg);
        if (!visible)
        {
            c.nextCheck = t + std::max(0.0, (m_minElevDeg - elev) / kMaxElevRate - 0.2);
        }
        if (!visible && !c.up)
        {
            continue;
        }
        const Vector p = c.mob->GetEcefPosition();
        const double dx = p.x - c.gsEcef.x;
        const double dy = p.y - c.gsEcef.y;
        const double dz = p.z - c.gsEcef.z;
        const double range = std::sqrt(dx * dx + dy * dy + dz * dz);
        ContactEvent ev{t, c.sat, c.gs, false, visible, range, elev};
        if (visible != c.up)
        {
            c.up = visible;
            m_gslState[{c.sat, c.gs}] = visible;
            if (visible)
            {
                ++m_gslUp;
                m_contactUp(ev);
            }
            else
            {
                ++m_gslDown;
                m_contactDown(ev);
            }
        }
        else if (visible)
        {
            m_contactUpdate(ev);
        }
    }

    // GSL visibility (full scan when no candidate list was set).
    if (m_gslCand.empty())
    for (auto& gsIt : m_gss)
    {""")
c = rep(c, """    if (!m_islCandidates.empty())
    {
        for (const auto& [a, b] : m_islCandidates)""", """    if (!m_islCandidates.empty() && (m_tickN++ % m_islEvery != 0))
    {
        // P7: not this tick
    }
    else if (!m_islCandidates.empty())
    {
        for (const auto& [a, b] : m_islCandidates)""")
W("contact-graph-scheduler.cc", c)

# ---------- router
h = R("contact-graph-router.h")
h = rep(h, "    /// Adjacency map: node -> set of neighbours.", """    /// P7: flat adjacency for Dijkstra and for weight lookups, neighbors kept in ascending order so ties
    /// resolve as in the map-based version. m_adj and m_edges stay for the rest of the interface.
    std::vector<std::vector<std::pair<uint32_t, double>>> m_flat;
    mutable std::vector<double> m_dist;
    mutable std::vector<uint32_t> m_pred;
    mutable std::vector<uint32_t> m_stamp;
    mutable uint32_t m_stampNow{0};
    void FlatAdd(uint32_t a, uint32_t b, double w);
    void FlatErase(uint32_t a, uint32_t b);
    void FlatSet(uint32_t a, uint32_t b, double w);

    /// Adjacency map: node -> set of neighbours.""")
W("contact-graph-router.h", h)
c = R("contact-graph-router.cc")
c = rep(c, 'NS_LOG_COMPONENT_DEFINE("ContactGraphRouter");', '''#ifdef P7_STANDALONE
NS_LOG_COMPONENT_DEFINE("ContactGraphRouterP7");
#else
NS_LOG_COMPONENT_DEFINE("ContactGraphRouter");
#endif''')
c = rep(c, 'TypeId("ns3::ntncon::ContactGraphRouter")', '''TypeId(
#ifdef P7_STANDALONE
                                   "ns3::ntncon::ContactGraphRouterP7"
#else
                                   "ns3::ntncon::ContactGraphRouter"
#endif
                                   )''')
c = rep(c, """void
ContactGraphRouter::UpdateEdgeWeight(const ContactEvent& ev)
{""", """void
ContactGraphRouter::FlatAdd(uint32_t a, uint32_t b, double w)
{
    if (m_flat.size() <= a)
    {
        m_flat.resize(a + 1);
    }
    auto& v = m_flat[a];
    auto it = std::lower_bound(v.begin(), v.end(), std::make_pair(b, -1e300));
    if (it != v.end() && it->first == b)
    {
        it->second = w;
    }
    else
    {
        v.insert(it, {b, w});
    }
}

void
ContactGraphRouter::FlatErase(uint32_t a, uint32_t b)
{
    if (a >= m_flat.size())
    {
        return;
    }
    auto& v = m_flat[a];
    v.erase(std::remove_if(v.begin(), v.end(), [b](const auto& x) { return x.first == b; }), v.end());
}

void
ContactGraphRouter::FlatSet(uint32_t a, uint32_t b, double w)
{
    if (a >= m_flat.size())
    {
        return;
    }
    for (auto& x : m_flat[a])
    {
        if (x.first == b)
        {
            x.second = w;
            return;
        }
    }
}

void
ContactGraphRouter::UpdateEdgeWeight(const ContactEvent& ev)
{""")
c = rep(c, """    const auto key = CanonicalEdge(ev.node_a, ev.node_b);
    if (m_edges.count(key) != 0)
    {
        m_edgeWeights[key] = ev.range_m;
    }""", """    FlatSet(ev.node_a, ev.node_b, ev.range_m);
    FlatSet(ev.node_b, ev.node_a, ev.range_m);""")
c = rep(c, """        // Seed the weight from the contact-up range; OnContactUpdate keeps
        // it current for the rest of the pass.
        m_edgeWeights[key] = ev.range_m;""", """        // Seed the weight from the contact-up range; OnContactUpdate keeps
        // it current for the rest of the pass.
        m_edgeWeights[key] = ev.range_m;
        FlatAdd(ev.node_a, ev.node_b, ev.range_m);
        FlatAdd(ev.node_b, ev.node_a, ev.range_m);""")
c = rep(c, """        m_edgeWeights.erase(key);
    }""", """        m_edgeWeights.erase(key);
        FlatErase(ev.node_a, ev.node_b);
        FlatErase(ev.node_b, ev.node_a);
    }""")
c = rep(c, """    auto it = m_edgeWeights.find(CanonicalEdge(a, b));
    if (it == m_edgeWeights.end())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return it->second;""", """    if (a < m_flat.size())
    {
        for (const auto& x : m_flat[a])
        {
            if (x.first == b)
            {
                return x.second;
            }
        }
    }
    return std::numeric_limits<double>::quiet_NaN();""")
c = rep(c, "    return m_edges.count(CanonicalEdge(a, b)) != 0;", """    return !std::isnan(EdgeWeight(a, b));""")
# flat Dijkstra
start = c.index("ContactGraphRouter::WeightedPath\nContactGraphRouter::ShortestPathWeightedFiltered")
end = c.index("std::vector<uint32_t>\nContactGraphRouter::ShortestPath(uint32_t src, uint32_t dst) const")
c = c[:start] + """ContactGraphRouter::WeightedPath
ContactGraphRouter::ShortestPathWeightedFiltered(uint32_t src,
                                                 uint32_t dst,
                                                 const EdgeFilter& keep) const
{
    ++m_queries;
    if (src == dst)
    {
        return {{src}, 0.0};
    }
    // Dijkstra over the flat adjacency; same relaxation order and tie rule as the map version.
    const size_t n = m_flat.size();
    if (src >= n || dst >= n)
    {
        return {{}, std::numeric_limits<double>::infinity()};
    }
    if (m_dist.size() < n)
    {
        m_dist.resize(n);
        m_pred.resize(n);
        m_stamp.assign(n, 0);
        m_stampNow = 0;
    }
    ++m_stampNow;
    const uint32_t st = m_stampNow;
    using QEntry = std::pair<double, uint32_t>; // (weight, node)
    std::priority_queue<QEntry, std::vector<QEntry>, std::greater<QEntry>> pq;
    m_dist[src] = 0.0;
    m_stamp[src] = st;
    pq.push({0.0, src});
    while (!pq.empty())
    {
        const auto [w, u] = pq.top();
        pq.pop();
        if (u == dst)
        {
            std::vector<uint32_t> path{dst};
            uint32_t c = dst;
            while (c != src)
            {
                c = m_pred[c];
                path.push_back(c);
            }
            std::reverse(path.begin(), path.end());
            return {std::move(path), w};
        }
        if (w > m_dist[u])
        {
            continue; // stale
        }
        for (const auto& [v, we] : m_flat[u])
        {
            if (keep && !keep(u, v))
            {
                continue;
            }
            const double newW = w + we;
            if (m_stamp[v] != st || newW < m_dist[v])
            {
                m_stamp[v] = st;
                m_dist[v] = newW;
                m_pred[v] = u;
                pq.push({newW, v});
            }
        }
    }
    return {{}, std::numeric_limits<double>::infinity()};
}

""" + c[end:]
c = rep(c, "#include <limits>\n", "#include <cmath>\n#include <limits>\n")
W("contact-graph-router.cc", c)
print("ok")
