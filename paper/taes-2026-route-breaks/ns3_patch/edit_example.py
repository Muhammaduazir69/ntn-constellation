import re
s=open('orig/ntn-constellation-route-lifetime.cc').read()
def rep(a,b,n=1):
    global s
    assert s.count(a)==n, (s.count(a),a)
    s=s.replace(a,b)
rep('#include "ns3/contact-graph-router.h"\n#include "ns3/contact-graph-scheduler.h"\n','#include "contact-graph-router.h"\n#include "contact-graph-scheduler.h"\n')
rep('#include <algorithm>\n','#include <sys/resource.h>\n#include <time.h>\n\n#include <algorithm>\n#include <fstream>\n')
rep(''' *   life      longest-lived path within (1 + --eps) of the shortest delay,
 *             switched --guard seconds before its predicted expiry
''',''' *   life      longest-lived path within (1 + --eps) of the shortest delay,
 *             switched --guard seconds before its predicted expiry
 * Added for paper 7 (Nov 2026):
 *   asrk      merged snapshot that keeps a route that still lasts the window
 *   track     Algorithm 1 on residual times estimated from a 30 s elevation track
 *   alpr      average-latency persistent routing of Bhattacharjee et al. (eq. 13)
 *   isasr     stability and activeness slotted routing of Bhattacharjee et al. (eq. 14-15)
''')
rep('''    uint64_t rx{0};
} g;
''','''    uint64_t rx{0};
    // ---- paper 7 additions
    bool noData{false};
    bool bbm{false};                   // one terminal per gateway: the old ground link stops when acquisition starts
    double fallbackEps{0.10};          // asr, asrk: stretch bound of the fallback route
    double keepEps{1e9};               // asrk: keep a route only while its stretch is within this
    double trackArc{30.0};             // track: arc of the elevation track (s)
    double gammaIsasr{1.0};            // isasr: gamma of eq. 14 (dimensionless)
    bool proactive{true};              // alpr, isasr: start the switch g + tAcq before the known end
    bool logEvents{false};
    uint64_t ingressChanges{0}, egressChanges{0}, islOnlyChanges{0}, acquisitions{0};
    uint64_t breaksIn{0}, breaksOut{0}, breaksIsl{0};
    uint64_t holdSwitches{0};
    double ageSum{0.0}, lastInstall{-1.0};
    uint64_t ageN{0};
    double resInSum{0.0}, resOutSum{0.0};
    uint64_t resN{0};
    double searchCpu{0.0};
    uint64_t searchN{0};
    std::set<Edge> blocked;
    std::map<Edge, double> trackCache;
    uint64_t trackEvals{0};
    double trackErrAbsSum{0.0};
    uint64_t trackErrN{0};
    double dropDeg{23.0};
} g;

double
ThreadCpu()
{
    timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

void Log(const std::string& what)
{
    if (g.logEvents)
    {
        std::cerr << Simulator::Now().GetSeconds() << " " << what << "\\n";
    }
}

std::string PathStr(const std::vector<uint32_t>& p)
{
    std::ostringstream o;
    for (auto v : p) o << v << " ";
    return o.str();
}
''')
rep('''void
InstallPath(const std::vector<uint32_t>& p)
{
    Ipv4StaticRoutingHelper srh;
''','''void Unblock();

void
CountChange(const std::vector<uint32_t>& p)
{
    const double now = Simulator::Now().GetSeconds();
    if (!g.path.empty() && p != g.path)
    {
        const bool in = p[1] != g.path[1];
        const bool out = p[p.size() - 2] != g.path[g.path.size() - 2];
        g.ingressChanges += in;
        g.egressChanges += out;
        g.islOnlyChanges += (!in && !out);
        g.acquisitions += (in || out);
        if (g.lastInstall >= 0)
        {
            g.ageSum += now - g.lastInstall;
            ++g.ageN;
        }
    }
    g.lastInstall = now;
    if (p.size() >= 3)
    {
        const double h = g.horizon;
        g.resInSum += std::min(h, g.predSched->PredictContactEndS(p[1], g.gsA, false, h, 5.0));
        g.resOutSum += std::min(h, g.predSched->PredictContactEndS(p[p.size() - 2], g.gsB, false, h, 5.0));
        ++g.resN;
    }
}

void
InstallPath(const std::vector<uint32_t>& p)
{
    Log("install " + PathStr(p));
    CountChange(p);
    if (g.noData)
    {
        if (!g.path.empty() && p != g.path)
        {
            ++g.changes;
        }
        g.path = p;
        g.pathEdges.clear();
        for (size_t k = 0; k + 1 < p.size(); ++k)
        {
            g.pathEdges.insert(Canon(p[k], p[k + 1]));
        }
        return;
    }
    Ipv4StaticRoutingHelper srh;
''')
rep('''        g.pathEdges.insert(Canon(p[k], p[k + 1]));
    }
}

bool
PathUp(''','''        g.pathEdges.insert(Canon(p[k], p[k + 1]));
    }
    Unblock();
}

bool
PathUp(''')
rep('''    g.pendingInstall.Cancel();
    if (g.tAcq > 0 && newGround && !g.path.empty())
    {
        g.pendingPath = p;
''','''    g.pendingInstall.Cancel();
    if (g.tAcq > 0 && newGround && !g.path.empty())
    {
        if (g.bbm && !g.noData)
        {
            // One terminal per gateway: the terminal that must slew stops serving its old link now.
            const Edge oldIn = Canon(g.path[1], g.gsA);
            const Edge oldOut = Canon(g.path[g.path.size() - 2], g.gsB);
            for (const Edge& oe : {oldIn, oldOut})
            {
                const Edge ne = oe.second == g.gsA ? Canon(p[1], g.gsA) : Canon(p[p.size() - 2], g.gsB);
                if (oe != ne && g.links.count(oe))
                {
                    g.links.at(oe).emA->SetRate(1.0);
                    g.links.at(oe).emB->SetRate(1.0);
                    g.blocked.insert(oe);
                }
            }
        }
        g.pendingPath = p;
''')
rep('''/// Elevation (deg) of satellite s seen from ground station gs (WGS84, geodetic up).
double
ElevDeg(uint32_t s, uint32_t gs)
{
    const int k = gs == g.gsA ? 0 : 1;''','''/// Elevation (deg) of satellite s seen from ground station gs (WGS84, geodetic up), offset s from now.
double
ElevDeg(uint32_t s, uint32_t gs, double offset = 0.0)
{
    const int k = gs == g.gsA ? 0 : 1;''')
rep('''    const Vector sp = g.satMob[s]->GetEcefPosition();
    const Vector d(sp.x - gp.x, sp.y - gp.y, sp.z - gp.z);''','''    Vector sp, sv;
    if (offset == 0.0)
    {
        sp = g.satMob[s]->GetEcefPosition();
    }
    else
    {
        g.satMob[s]->GetEcefStateAt(Seconds(offset), sp, sv);
    }
    const Vector d(sp.x - gp.x, sp.y - gp.y, sp.z - gp.z);''')
rep('''/// Elevation-trigger policy (no ephemeris): every dt,''','''/// Central angle (rad) from gateway to sub-satellite point at elevation e (deg), for shell radius r.
double
CentralAngle(double eDeg, double r)
{
    const double e = eDeg * M_PI / 180;
    return std::acos(6371000.0 * std::cos(e) / r) - e;
}

/// Residual time of a ground link estimated from three elevation samples over the last trackArc seconds, with no
/// orbit: psi^2 of a pass is close to a parabola in time (the satellite crosses the visibility cap along a chord at
/// constant speed), so a parabola through psi^2 at the three samples gives the time at which psi reaches the drop
/// elevation. Inter-satellite links are taken not to end (the tracker cannot see them).
double
TrackResidual(uint32_t a, uint32_t b)
{
    const Edge e = Canon(a, b);
    if (e.second < g.n)
    {
        return g.horizon;
    }
    auto it = g.trackCache.find(e);
    if (it != g.trackCache.end())
    {
        return it->second;
    }
    ++g.trackEvals;
    const uint32_t s = e.first, gs = e.second;
    const double h = g.trackArc / 2;
    Vector sp = g.satMob[s]->GetEcefPosition();
    const double r = std::sqrt(sp.x * sp.x + sp.y * sp.y + sp.z * sp.z);
    const double p0 = CentralAngle(ElevDeg(s, gs), r);
    const double p1 = CentralAngle(ElevDeg(s, gs, -h), r);
    const double p2 = CentralAngle(ElevDeg(s, gs, -2 * h), r);
    const double y0 = p0 * p0, y1 = p1 * p1, y2 = p2 * p2;
    const double C = (y2 - 2 * y1 + y0) / (2 * h * h);
    const double B = (y0 - y1) / h + C * h;
    const double pd = CentralAngle(g.dropDeg, r);
    const double A = y0 - pd * pd; // A + B t + C t^2 = 0 at the drop
    double res = 0.0;
    if (A < 0 && C > 0)
    {
        res = (-B + std::sqrt(std::max(0.0, B * B - 4 * A * C))) / (2 * C);
    }
    else if (A < 0)
    {
        res = B > 0 ? -A / B : g.horizon;
    }
    res = std::min(res, g.horizon);
    g.trackCache[e] = res;
    return res;
}

/// One-way delay (m of path length) of route p at the given offset from now, from the ephemeris.
double
RouteLengthAt(const std::vector<uint32_t>& p, double offset)
{
    auto posOf = [&](uint32_t v) {
        Vector x, u;
        if (v < g.n)
        {
            if (offset == 0.0)
            {
                x = g.satMob[v]->GetEcefPosition();
            }
            else
            {
                g.satMob[v]->GetEcefStateAt(Seconds(offset), x, u);
            }
            return x;
        }
        const int k = v == g.gsA ? 0 : 1;
        const double la = g.gsLat[k] * M_PI / 180, lo = g.gsLon[k] * M_PI / 180;
        const double a = 6378137.0, e2 = 6.69437999014e-3;
        const double nn = a / std::sqrt(1 - e2 * std::sin(la) * std::sin(la));
        return Vector(nn * std::cos(la) * std::cos(lo), nn * std::cos(la) * std::sin(lo), nn * (1 - e2) * std::sin(la));
    };
    double L = 0;
    Vector prev = posOf(p[0]);
    for (size_t k = 1; k < p.size(); ++k)
    {
        const Vector c = posOf(p[k]);
        const double dx = c.x - prev.x, dy = c.y - prev.y, dz = c.z - prev.z;
        L += std::sqrt(dx * dx + dy * dy + dz * dz);
        prev = c;
    }
    return L;
}

/// Elevation-trigger policy (no ephemeris): every dt,''')
rep('''    ++g.computations;
    std::vector<uint32_t> p;
    if (g.policy == "life")
    {
        auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &Residual, g.guard + g.tAcq);
        g.dijkstra += lp.dijkstra_runs;
''','''    ++g.computations;
    std::vector<uint32_t> p;
    auto lifeOf = [](const std::vector<uint32_t>& q) {
        double l = g.horizon;
        for (size_t k = 0; k + 1 < q.size(); ++k)
        {
            l = std::min(l, Residual(q[k], q[k + 1]));
        }
        return l;
    };
    if (g.policy == "track")
    {
        // Algorithm 1 fed by the elevation track, switched by the elevation trigger (ElevCheck)
        g.trackCache.clear();
        const double t0 = ThreadCpu();
        auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &TrackResidual, g.guard + g.tAcq);
        g.searchCpu += ThreadCpu() - t0;
        ++g.searchN;
        g.dijkstra += lp.dijkstra_runs;
        p = lp.path;
    }
    else if (g.policy == "alpr")
    {
        // Bhattacharjee et al., Algorithm 2 and eq. (13) with eta_s = tAcq and 1 s slots: edge-disjoint shortest routes,
        // each scored by (eta_s + sum of delays over the slots it exists) / slots; the lowest is kept until it ends.
        std::set<Edge> removed;
        std::vector<std::vector<uint32_t>> cand;
        const size_t nr = std::min(g.router->Neighbours(g.gsA).size(), g.router->Neighbours(g.gsB).size());
        for (size_t k = 0; k < nr; ++k)
        {
            auto wp = g.router->ShortestPathWeightedFiltered(g.gsA, g.gsB, [&](uint32_t a, uint32_t b) {
                return !removed.count(Canon(a, b));
            });
            ++g.dijkstra;
            if (wp.path.empty())
            {
                break;
            }
            for (size_t q = 0; q + 1 < wp.path.size(); ++q)
            {
                removed.insert(Canon(wp.path[q], wp.path[q + 1]));
            }
            cand.push_back(wp.path);
        }
        double best = std::numeric_limits<double>::infinity();
        double bestL = 0;
        for (const auto& c : cand)
        {
            const double l = lifeOf(c);
            if (l <= g.guard + g.tAcq && cand.size() > 1 && c != cand.back())
            {
                continue; // cannot be switched to in time
            }
            double sum = 0;
            for (int q = 0; q < 5; ++q)
            {
                sum += RouteLengthAt(c, l * (q + 0.5) / 5) / kC;
            }
            const double score = (g.tAcq + (sum / 5) * std::max(1.0, std::floor(l))) / std::max(1.0, std::floor(l));
            if (score < best)
            {
                best = score;
                p = c;
                bestL = l;
            }
        }
        if (!p.empty() && g.proactive)
        {
            const double wait = bestL - g.guard - g.tAcq;
            Schedule(wait > 1e-3 ? wait : 0.5);
        }
    }
    else if (g.policy == "isasr")
    {
        // Bhattacharjee et al., Algorithm 3 and eq. (14)-(15) with eta_s = tAcq, 1 s slots, gamma = gammaIsasr, on the
        // +Grid topology: an edge of the active route has no activeness cost; an inter-satellite link that outlives the
        // horizon has no stability cost and, being permanent, no activeness cost either. Until the active route's break
        // point its edges keep cost_act = 0; at the break point (here g + tAcq before it) all edges carry eta_s.
        const double etaS = std::max(g.tAcq, 1e-3);
        const double lActive = PathUp(g.path) ? lifeOf(g.path) : 0.0;
        const bool atBreak = !g.path.empty() && lActive <= g.guard + g.tAcq + 1.0; // one slot of margin
        const std::set<Edge> active = atBreak ? std::set<Edge>{} : g.pathEdges;
        auto wp = g.router->ShortestPathCost(g.gsA, g.gsB, [&](uint32_t a, uint32_t b, double w) {
            const Edge e = Canon(a, b);
            const double res = Residual(a, b);
            if (e.second < g.n && res >= g.horizon)
            {
                return w;
            }
            if (res <= g.guard + g.tAcq && !active.count(e))
            {
                return std::numeric_limits<double>::infinity();
            }
            const double costSt = etaS / std::max(1.0, std::floor(res));
            const double costAct = active.count(e) ? 0.0 : etaS;
            return w + g.gammaIsasr * kC * (costSt + costAct); // weights are in meters of path; costs converted by c
        });
        ++g.dijkstra;
        p = wp.path;
        Schedule(1.0);
    }
    else if (g.policy == "life")
    {
        const double t0 = ThreadCpu();
        auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &Residual, g.guard + g.tAcq);
        g.searchCpu += ThreadCpu() - t0;
        ++g.searchN;
        g.dijkstra += lp.dijkstra_runs;
''')
rep('''    else if (g.policy == "asr")
    {
        auto wp = g.router->ShortestPathWeightedFiltered''','''    else if (g.policy == "asr" || g.policy == "asrk")
    {
        if (g.policy == "asrk" && PathUp(g.path) && lifeOf(g.path) >= g.window + g.tAcq)
        {
            // keep a route that still lasts the window, while its stretch is within the keep bound
            double w = 0;
            for (size_t k = 0; k + 1 < g.path.size(); ++k)
            {
                w += g.router->EdgeWeight(g.path[k], g.path[k + 1]);
            }
            const double best = g.router->ShortestPathWeighted(g.gsA, g.gsB).total_weight;
            ++g.dijkstra;
            if (w <= (1.0 + g.keepEps) * best)
            {
                Schedule(g.window);
                return;
            }
        }
        auto wp = g.router->ShortestPathWeightedFiltered''')
rep('''            // nothing lasts the window: fall back to the longest-lived route, and repair it on a break (OnDown)
            auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, 1e9, &Residual, g.guard + g.tAcq);''','''            // nothing lasts the window: fall back to the longest-lived route within fallbackEps, repair on a break (OnDown)
            auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.fallbackEps, &Residual, g.guard + g.tAcq);''')
rep('''    if (g.path.empty())
    {
        InstallPath(p); // the first route is in place before traffic starts''','''    Log("decide " + PathStr(p));
    if (g.path.empty())
    {
        InstallPath(p); // the first route is in place before traffic starts''')
rep('''    if (!g.pathEdges.count(e))
    {
        return;
    }
    ++g.pathBreaks;
''','''    if (!g.pathEdges.count(e))
    {
        return;
    }
    ++g.pathBreaks;
    (e.second == g.gsA ? g.breaksIn : e.second == g.gsB ? g.breaksOut : g.breaksIsl)++;
    Log("break " + std::to_string(e.first) + "-" + std::to_string(e.second));
''')
rep('''    if (g.policy == "ondemand" || g.policy == "elev" || g.policy == "elevb" || g.policy == "asr")
''','''    if (g.policy == "ondemand" || g.policy == "elev" || g.policy == "elevb" || g.policy == "asr" ||
        g.policy == "asrk" || g.policy == "track" || g.policy == "alpr")
''')
rep('''    if (it != g.links.end())
    {
        it->second.emA->SetRate(0.0);
        it->second.emB->SetRate(0.0);
    }
}
''','''    if (it != g.links.end() && !g.blocked.count(e))
    {
        it->second.emA->SetRate(0.0);
        it->second.emB->SetRate(0.0);
    }
}

/// The acquisition has finished: terminals that were blocked serve their links again if those are still up.
void
Unblock()
{
    for (const Edge& e : g.blocked)
    {
        auto it = g.links.find(e);
        if (it != g.links.end() && g.router->HasEdge(e.first, e.second))
        {
            it->second.emA->SetRate(0.0);
            it->second.emB->SetRate(0.0);
        }
    }
    g.blocked.clear();
}
''')
rep('''            if (g.hold)
            {
                ++g.dijkstra; // the hold check needs this shortest path
                if (w / best - 1.0 > g.eps)
                {
                    Schedule(0.0);
                }
            }''','''            if (g.hold)
            {
                ++g.dijkstra; // the hold check needs this shortest path
                if (w / best - 1.0 > g.eps && !g.pendingInstall.IsPending())
                {
                    ++g.holdSwitches;
                    Schedule(0.0);
                }
            }''')
rep('''    if (g.policy == "elev" || g.policy == "elevb")
    {''','''    if (g.policy == "elev" || g.policy == "elevb" || g.policy == "track")
    {''')
rep('''    cmd.AddValue("reachOut", "write the longest admissible lifetime every second to this CSV", reachOut);
''','''    cmd.AddValue("reachOut", "write the longest admissible lifetime every second to this CSV", reachOut);
    cmd.AddValue("noData", "control plane only: no IP stack, no packets (same decisions, loss not measured)", g.noData);
    cmd.AddValue("bbm", "one terminal per gateway: old ground link stops when acquisition starts", g.bbm);
    cmd.AddValue("fallbackEps", "asr, asrk: stretch bound of the longest-lived fallback route", g.fallbackEps);
    cmd.AddValue("keepEps", "asrk: keep a route only while its stretch is within this", g.keepEps);
    cmd.AddValue("trackArc", "track: arc of the elevation track (s)", g.trackArc);
    cmd.AddValue("gamma", "isasr: gamma of eq. (14), dimensionless", g.gammaIsasr);
    cmd.AddValue("proactive", "alpr: start the switch before the known end of the route", g.proactive);
    cmd.AddValue("logEvents", "print decisions, installs and breaks to stderr", g.logEvents);
    double islRefresh = 1.0;
    cmd.AddValue("islRefresh", "refresh interval of inter-satellite ranges (s), 0 = every tick", islRefresh);
    std::string tleFile, edgeFile;
    double epochUnix = 1792800000.0;
    cmd.AddValue("tle", "real shell: TLE file (3-line sets), propagated with Vallado SGP4", tleFile);
    cmd.AddValue("edges", "real shell: +Grid candidates, lines 'a b cross'", edgeFile);
    cmd.AddValue("epochUnix", "Unix time of simulation start before startOffset", epochUnix);
''')
rep('''    wc.epoch_unix_s = 1792800000.0; // 2026-10-24
    const auto els = star ? WalkerConstellation::BuildStar(wc) : WalkerConstellation::BuildDelta(wc);
    g.n = els.size();''','''    wc.epoch_unix_s = 1792800000.0; // 2026-10-24
    auto els = star ? WalkerConstellation::BuildStar(wc) : WalkerConstellation::BuildDelta(wc);
    std::vector<ns3::ntncon::TleRecord> tles;
    if (!tleFile.empty())
    {
        std::ifstream tf(tleFile);
        std::stringstream body;
        body << tf.rdbuf();
        ns3::ntncon::ParseTleStream(body.str(), tles);
        els.assign(tles.size(), els[0]);
    }
    g.n = els.size();''')
rep('''        auto m = CreateObject<Sgp4MobilityModel>();
        m->SetElements(els[k]);
        m->SetSimulationEpochUnix(wc.epoch_unix_s + startOffset);
        sats.push_back(m);''','''        auto m = CreateObject<Sgp4MobilityModel>();
        if (tles.empty())
        {
            m->SetElements(els[k]);
            m->SetSimulationEpochUnix(wc.epoch_unix_s + startOffset);
        }
        else
        {
            m->SetTle(tles[k]);
            m->SetUseVallado(true);
            NS_ABORT_MSG_IF(!m->IsUsingSgp4(), "SGP4 init failed for satellite " << k);
            m->SetSimulationEpochUnix(epochUnix + startOffset);
        }
        sats.push_back(m);''')
rep('''    for (uint32_t p = 0; p < planes; ++p)
    {
        for (uint32_t s = 0; s < perPlane; ++s)
        {
            const uint32_t a = p * perPlane + s;''','''    if (!edgeFile.empty())
    {
        std::ifstream ef(edgeFile);
        uint32_t a, b;
        int isCross;
        while (ef >> a >> b >> isCross)
        {
            isl.push_back(Canon(a, b));
            if (isCross)
            {
                cross.insert(isl.back());
            }
        }
    }
    else
    for (uint32_t p = 0; p < planes; ++p)
    {
        for (uint32_t s = 0; s < perPlane; ++s)
        {
            const uint32_t a = p * perPlane + s;''')
rep('''    // ---- data plane
    g.nodes.Create(g.n + 2);
    InternetStackHelper stack;
    Ipv4StaticRoutingHelper srh;
    stack.SetRoutingHelper(srh);
    stack.Install(g.nodes);
''','''    g.sched->SetGslCandidates(gsl);
    // Inter-satellite ranges move slowly and no +Grid link drops under the range cap, so they are refreshed every
    // islRefresh seconds; with a latitude cutoff the links do drop, so those runs sample on every tick.
    g.sched->SetIslRefreshEvery(crossLat > 0 || islRefresh <= 0 ? 1 : static_cast<uint32_t>(std::lround(islRefresh / tick)));
    g.dropDeg = minElev - 2.0;
    g.horizon = 1800.0;

    // ---- data plane
    InternetStackHelper stack;
    Ipv4StaticRoutingHelper srh;
    if (!g.noData)
    {
    g.nodes.Create(g.n + 2);
    stack.SetRoutingHelper(srh);
    stack.Install(g.nodes);
    }
''')
rep('''    for (const auto& e : isl)
    {
        mk(e);
    }
    for (const auto& e : gsl)
    {
        mk(e);
    }
''','''    if (!g.noData)
    {
    for (const auto& e : isl)
    {
        mk(e);
    }
    for (const auto& e : gsl)
    {
        mk(e);
    }
    }
''')
rep('''    const uint16_t port = 9000;
    PacketSinkHelper sinkH("ns3::UdpSocketFactory", InetSocketAddress(Ipv4Address::GetAny(), port));
    auto sinkApp = sinkH.Install(g.nodes.Get(g.gsB));
    sinkApp.Get(0)->TraceConnectWithoutContext("Rx", MakeCallback(&OnRx));
    sinkApp.Start(Seconds(0));
    UdpClientHelper cli(g.dst, port);
    const double t0 = 5.0;
    const uint64_t nPk = static_cast<uint64_t>((duration - t0) * rate);
    cli.SetAttribute("MaxPackets", UintegerValue(nPk));
    cli.SetAttribute("Interval", TimeValue(Seconds(1.0 / rate)));
    cli.SetAttribute("PacketSize", UintegerValue(200));
    auto cliApp = cli.Install(g.nodes.Get(g.gsA));
    cliApp.Start(Seconds(t0));
    cliApp.Stop(Seconds(duration));
''','''    const uint16_t port = 9000;
    const double t0 = 5.0;
    const uint64_t nPk = static_cast<uint64_t>((duration - t0) * rate);
    if (!g.noData)
    {
    PacketSinkHelper sinkH("ns3::UdpSocketFactory", InetSocketAddress(Ipv4Address::GetAny(), port));
    auto sinkApp = sinkH.Install(g.nodes.Get(g.gsB));
    sinkApp.Get(0)->TraceConnectWithoutContext("Rx", MakeCallback(&OnRx));
    sinkApp.Start(Seconds(0));
    UdpClientHelper cli(g.dst, port);
    cli.SetAttribute("MaxPackets", UintegerValue(nPk));
    cli.SetAttribute("Interval", TimeValue(Seconds(1.0 / rate)));
    cli.SetAttribute("PacketSize", UintegerValue(200));
    auto cliApp = cli.Install(g.nodes.Get(g.gsA));
    cliApp.Start(Seconds(t0));
    cliApp.Stop(Seconds(duration));
    }
''')
rep('''        g.changes = g.pathBreaks = g.computations = g.dijkstra = g.predictions = g.forced = 0;
''','''        g.changes = g.pathBreaks = g.computations = g.dijkstra = g.predictions = g.forced = 0;
        g.ingressChanges = g.egressChanges = g.islOnlyChanges = g.acquisitions = 0;
        g.breaksIn = g.breaksOut = g.breaksIsl = g.holdSwitches = g.ageN = g.resN = g.searchN = 0;
        g.ageSum = g.resInSum = g.resOutSum = g.searchCpu = 0.0;
''')
rep('''    js << "{\\"ver\\":3,\\"policy\\":\\"" << g.policy''','''    rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    auto shareOver = [](double x) {
        return g.stretch.empty() ? NAN
                                 : static_cast<double>(std::count_if(g.stretch.begin(), g.stretch.end(),
                                                                     [x](double v) { return v > x; })) /
                                       g.stretch.size();
    };
    js << "{\\"ver\\":7,\\"bbm\\":" << g.bbm << ",\\"no_data\\":" << g.noData
       << ",\\"ingress_per_h\\":" << g.ingressChanges / hours << ",\\"egress_per_h\\":" << g.egressChanges / hours
       << ",\\"isl_only_per_h\\":" << g.islOnlyChanges / hours << ",\\"acq_per_h\\":" << g.acquisitions / hours
       << ",\\"breaks_in_per_h\\":" << g.breaksIn / hours << ",\\"breaks_out_per_h\\":" << g.breaksOut / hours
       << ",\\"breaks_isl_per_h\\":" << g.breaksIsl / hours << ",\\"hold_per_h\\":" << g.holdSwitches / hours
       << ",\\"age_mean_s\\":" << (g.ageN ? g.ageSum / g.ageN : NAN)
       << ",\\"res_in_mean_s\\":" << (g.resN ? g.resInSum / g.resN : NAN)
       << ",\\"res_out_mean_s\\":" << (g.resN ? g.resOutSum / g.resN : NAN)
       << ",\\"search_cpu_ms\\":" << (g.searchN ? 1e3 * g.searchCpu / g.searchN : NAN)
       << ",\\"stretch_p50\\":" << Pct(g.stretch, 0.5) << ",\\"stretch_p95\\":" << Pct(g.stretch, 0.95)
       << ",\\"stretch_max\\":" << Pct(g.stretch, 1.0) << ",\\"over5\\":" << shareOver(0.05)
       << ",\\"over10\\":" << shareOver(0.10) << ",\\"cpu_s\\":" << ru.ru_utime.tv_sec + 1e-6 * ru.ru_utime.tv_usec
       << ",\\"propagator\\":\\"" << (sats[0]->IsUsingSgp4() ? "sgp4" : "kepler-j2") << "\\"";
    js << ",\\"policy\\":\\"" << g.policy''')

# ---- admissibility: a route must outlive its own acquisition, its successor's acquisition and the guard (M9)
rep("""    auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &Residual, g.guard + g.tAcq);
    const double best = g.router->ShortestPathWeighted(g.gsA, g.gsB).total_weight;""","""    auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &Residual, g.guard + 2 * g.tAcq);
    const double best = g.router->ShortestPathWeighted(g.gsA, g.gsB).total_weight;""")
rep("""auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &TrackResidual, g.guard + g.tAcq);""","""auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &TrackResidual, g.guard + 2 * g.tAcq);""")
rep("""        auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &Residual, g.guard + g.tAcq);
        g.searchCpu += ThreadCpu() - t0;
        ++g.searchN;
        g.dijkstra += lp.dijkstra_runs;
""","""        auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.eps, &Residual, g.guard + 2 * g.tAcq);
        g.searchCpu += ThreadCpu() - t0;
        ++g.searchN;
        g.dijkstra += lp.dijkstra_runs;
        if (lp.lifetime_s <= g.guard + 2 * g.tAcq)
        {
            ++g.noAdmissible;
        }
""")
rep("""auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.fallbackEps, &Residual, g.guard + g.tAcq);""","""auto lp = g.router->LifetimeAwarePath(g.gsA, g.gsB, g.fallbackEps, &Residual, g.guard + 2 * g.tAcq);""")
rep("""            if (l <= g.guard + g.tAcq && cand.size() > 1 && c != cand.back())""","""            if (l <= g.guard + 2 * g.tAcq && cand.size() > 1 && c != cand.back())""")
rep("""            if (res <= g.guard + g.tAcq && !active.count(e))""","""            if (res <= g.guard + 2 * g.tAcq && !active.count(e))""")
rep("""    uint64_t trackEvals{0};
""","""    uint64_t trackEvals{0};
    uint64_t noAdmissible{0};
    double routeLenSum{0.0};
    uint64_t routeLenN{0};
""")
rep("""            g.stretch.push_back(w / best - 1.0);
""","""            g.stretch.push_back(w / best - 1.0);
            g.routeLenSum += w;
            ++g.routeLenN;
""")
rep("""       << ",\\"propagator\\":\\"" <<""","""       << ",\\"no_admissible\\":" << g.noAdmissible
       << ",\\"route_ms_mean\\":" << (g.routeLenN ? 1e3 * g.routeLenSum / g.routeLenN / kC : NAN)
       << ",\\"propagator\\":\\"" <<""")
rep("""        g.breaksIn = g.breaksOut = g.breaksIsl = g.holdSwitches = g.ageN = g.resN = g.searchN = 0;
""","""        g.breaksIn = g.breaksOut = g.breaksIsl = g.holdSwitches = g.ageN = g.resN = g.searchN = g.noAdmissible = 0;
""")

# ==== round 2 of review ====
# C1: a link that outlives the prediction horizon must not age toward a false end. Saturated predictions are
# remembered with the time they were made and re-evaluated after satRecheck seconds.
rep("""    auto it = g.endCache.find(e);
    if (it != g.endCache.end() && it->second > now)
    {
        return it->second - now;
    }
    const bool isl = e.second < g.n;
    ++g.predictions;
    const double r =
        g.predSched->PredictContactEndS(e.first, e.second, isl, g.horizon, isl ? 30.0 : 5.0);
    g.endCache[e] = now + r;
    return r;""","""    auto sat = g.satCache.find(e);
    if (sat != g.satCache.end() && now - sat->second < g.satRecheck)
    {
        return g.horizon; // no end within the horizon at the last look; it does not age
    }
    auto it = g.endCache.find(e);
    if (it != g.endCache.end() && it->second > now)
    {
        return it->second - now;
    }
    const bool isl = e.second < g.n;
    ++g.predictions;
    const double r =
        g.predSched->PredictContactEndS(e.first, e.second, isl, g.horizon, isl ? 30.0 : 5.0);
    if (r >= g.horizon)
    {
        g.satCache[e] = now;
        g.endCache.erase(e);
    }
    else
    {
        g.satCache.erase(e);
        g.endCache[e] = now + r;
    }
    return r;""")
rep("""    std::map<Edge, double> endCache; // absolute predicted contact end (s)
""","""    std::map<Edge, double> endCache; // absolute predicted contact end (s)
    std::map<Edge, double> satCache; // link with no end within the horizon: time of that prediction
    double satRecheck{900.0};
""")
rep("""void
OnUp(ContactEvent ev)
{
    const Edge e = Canon(ev.node_a, ev.node_b);
    g.endCache.erase(e);
""","""void
OnUp(ContactEvent ev)
{
    const Edge e = Canon(ev.node_a, ev.node_b);
    g.endCache.erase(e);
    g.satCache.erase(e);
""")
rep("""void
OnDown(ContactEvent ev)
{
    const Edge e = Canon(ev.node_a, ev.node_b);
    g.endCache.erase(e);
""","""void
OnDown(ContactEvent ev)
{
    const Edge e = Canon(ev.node_a, ev.node_b);
    g.endCache.erase(e);
    g.satCache.erase(e);
""")
# C2: ISASR as Algorithm 3 of the source: cost_act is state per edge, zero after an edge served an active route until
# that route's break point, eta_s otherwise; only edges of the active route are updated. Links that do not end within the
# horizon are permanent +Grid links: no stability cost (l = N in eq. 15) and, being up already, no activeness cost.
rep("""        const double lActive = PathUp(g.path) ? lifeOf(g.path) : 0.0;
        const bool atBreak = !g.path.empty() && lActive <= g.guard + g.tAcq + 1.0; // one slot of margin
        const std::set<Edge> active = atBreak ? std::set<Edge>{} : g.pathEdges;
""","""        const double lActive = PathUp(g.path) ? lifeOf(g.path) : 0.0;
        const bool atBreak = !g.path.empty() && lActive <= g.guard + g.tAcq + 1.0; // one slot of margin
        if (atBreak)
        {
            for (const Edge& pe : g.pathEdges)
            {
                g.isasrZeroAct.erase(pe);
            }
        }
        const std::set<Edge>& active = g.isasrZeroAct;
""")
rep("""        ++g.dijkstra;
        p = wp.path;
        Schedule(1.0);
""","""        ++g.dijkstra;
        p = wp.path;
        if (!p.empty() && lifeOf(p) > g.guard + g.tAcq + 1.0)
        {
            for (size_t q = 0; q + 1 < p.size(); ++q)
            {
                g.isasrZeroAct.insert(Canon(p[q], p[q + 1]));
            }
        }
        Schedule(1.0);
""")
rep("""    uint64_t trackEvals{0};
""","""    uint64_t trackEvals{0};
    std::set<Edge> isasrZeroAct;
    double elevNoiseDeg{0.0};
    uint32_t trackSamples{3};
    double floorGrid{5.0};
    Ptr<NormalRandomVariable> noise;
    std::map<std::pair<uint32_t, uint32_t>, double> elevReading; // one reading per link per decision
    uint64_t trackBad{0};
    uint64_t fallbackUsed{0};
    uint64_t mergedDecisions{0};
""")
# elevation as measured: truth plus noise, one reading per link per decision
rep("""/// Central angle (rad) from gateway""","""/// Elevation as a gateway would measure it now.
double
MeasElev(uint32_t s, uint32_t gs, double offset = 0.0)
{
    const double t = ElevDeg(s, gs, offset);
    return g.elevNoiseDeg > 0 ? t + g.noise->GetValue(0.0, g.elevNoiseDeg * g.elevNoiseDeg) : t;
}

/// The reading of one link within one decision.
double
ReadElev(uint32_t s, uint32_t gs)
{
    if (g.elevNoiseDeg <= 0)
    {
        return ElevDeg(s, gs);
    }
    auto key = std::make_pair(s, gs);
    auto it = g.elevReading.find(key);
    if (it != g.elevReading.end())
    {
        return it->second;
    }
    const double v = MeasElev(s, gs);
    g.elevReading[key] = v;
    return v;
}

/// Central angle (rad) from gateway""")
a=s.index("    const double h = g.trackArc / 2;\n")
b=s.index("    res = std::min(res, g.horizon);\n    g.trackCache[e] = res;")
s=s[:a]+"""    // Least-squares parabola y(tau) = c0 + c1 tau + c2 tau^2 through trackSamples readings of psi^2 over the arc;
    // with three samples it passes through them.
    Vector sp = g.satMob[s]->GetEcefPosition();
    const double r = std::sqrt(sp.x * sp.x + sp.y * sp.y + sp.z * sp.z);
    const uint32_t N = std::max<uint32_t>(3, g.trackSamples);
    double S[5] = {0, 0, 0, 0, 0}, Y[3] = {0, 0, 0};
    for (uint32_t k = 0; k < N; ++k)
    {
        const double tau = -g.trackArc * k / (N - 1);
        const double p = CentralAngle(MeasElev(s, gs, tau), r);
        const double y = p * p;
        double tk = 1.0;
        for (int q = 0; q < 5; ++q)
        {
            S[q] += tk;
            if (q < 3)
            {
                Y[q] += tk * y;
            }
            tk *= tau;
        }
    }
    // normal equations, Cramer
    const double D = S[0] * (S[2] * S[4] - S[3] * S[3]) - S[1] * (S[1] * S[4] - S[3] * S[2]) + S[2] * (S[1] * S[3] - S[2] * S[2]);
    const double c0 = (Y[0] * (S[2] * S[4] - S[3] * S[3]) - S[1] * (Y[1] * S[4] - S[3] * Y[2]) + S[2] * (Y[1] * S[3] - S[2] * Y[2])) / D;
    const double B = (S[0] * (Y[1] * S[4] - S[3] * Y[2]) - Y[0] * (S[1] * S[4] - S[3] * S[2]) + S[2] * (S[1] * Y[2] - Y[1] * S[2])) / D;
    const double C = (S[0] * (S[2] * Y[2] - Y[1] * S[3]) - S[1] * (S[1] * Y[2] - Y[1] * S[2]) + Y[0] * (S[1] * S[3] - S[2] * S[2])) / D;
    const double pd = CentralAngle(g.dropDeg, r);
    const double A = c0 - pd * pd; // A + B t + C t^2 = 0 at the drop
    double res = 0.0;
    if (A < 0 && C > 0)
    {
        res = (-B + std::sqrt(std::max(0.0, B * B - 4 * A * C))) / (2 * C);
    }
    else if (A < 0)
    {
        ++g.trackBad; // the arc does not fix a pass: a falling reading extrapolates linearly, a rising one is taken as lasting
        res = B > 0 ? -A / B : g.horizon;
    }
"""+s[b:]
# trigger and rank use readings
rep("""        if (ElevDeg(sA, g.gsA) < g.trigger || ElevDeg(sB, g.gsB) < g.trigger)""","""        if (MeasElev(sA, g.gsA) < g.trigger || MeasElev(sB, g.gsB) < g.trigger)""")
rep("""            return ElevDeg(std::min(a, b), std::max(a, b)) >= g.trigger;
        });
        ++g.dijkstra;
        p = base.path;""","""            return ReadElev(std::min(a, b), std::max(a, b)) >= g.trigger;
        });
        ++g.dijkstra;
        p = base.path;""")
rep("""            for (double floorDeg = 85.0; floorDeg > g.trigger; floorDeg -= 5.0)""","""            for (double floorDeg = 85.0; floorDeg > g.trigger; floorDeg -= g.floorGrid)""")
rep("""                    return ElevDeg(std::min(a, b), std::max(a, b)) >= floorDeg;""","""                    return ReadElev(std::min(a, b), std::max(a, b)) >= floorDeg;""")
rep("""            return ElevDeg(s, gs) >= g.trigger;
        });""","""            return ReadElev(s, gs) >= g.trigger;
        });""")
rep("""    ++g.computations;
    std::vector<uint32_t> p;
    auto lifeOf""","""    ++g.computations;
    g.elevReading.clear();
    std::vector<uint32_t> p;
    auto lifeOf""")
# merged fallback share (m6)
rep("""        ++g.dijkstra;
        if (wp.path.empty())
        {
            // nothing lasts the window""","""        ++g.dijkstra;
        ++g.mergedDecisions;
        if (wp.path.empty())
        {
            ++g.fallbackUsed;
            // nothing lasts the window""")
rep("""    cmd.AddValue("logEvents",""","""    cmd.AddValue("elevNoiseDeg", "standard deviation of measured elevation per reading (deg)", g.elevNoiseDeg);
    cmd.AddValue("trackSamples", "track: readings over the arc, fitted by least squares", g.trackSamples);
    cmd.AddValue("floorGrid", "elevb: grid of the elevation floor (deg)", g.floorGrid);
    cmd.AddValue("logEvents",""")
rep("""    g.sched->SetGslCandidates(gsl);
""","""    g.sched->SetGslCandidates(gsl);
    g.noise = CreateObject<NormalRandomVariable>();
""")
rep("""       << ",\\"no_admissible\\":" << g.noAdmissible
""","""       << ",\\"no_admissible\\":" << g.noAdmissible << ",\\"track_bad\\":" << g.trackBad
       << ",\\"fallback_share\\":" << (g.mergedDecisions ? static_cast<double>(g.fallbackUsed) / g.mergedDecisions : NAN)
""")
rep("""g.searchN = g.noAdmissible = 0;""","""g.searchN = g.noAdmissible = g.trackBad = g.fallbackUsed = g.mergedDecisions = 0;""")

# ==== round 3 of review: ISASR activeness on every edge, as in Algorithm 3 of the source (flag; default off keeps r2 behavior)
rep("""            if (e.second < g.n && res >= g.horizon)
            {
                return w;
            }
            if (res <= g.guard + 2 * g.tAcq && !active.count(e))""","""            if (e.second < g.n && res >= g.horizon)
            {
                // permanent +Grid link: l = N in eq. (15), so no stability cost. With isasrIslAct its activeness cost follows
                // Algorithm 3 (eta_s until it serves an active route); without it, as in our variant, it costs nothing.
                return g.isasrIslAct && !active.count(e) ? w + g.gammaIsasr * kC * etaS : w;
            }
            if (res <= g.guard + 2 * g.tAcq && !active.count(e))""")
rep("""    std::set<Edge> isasrZeroAct;
""","""    std::set<Edge> isasrZeroAct;
    bool isasrIslAct{false};
""")
rep("""    cmd.AddValue("gamma", "isasr: gamma of eq. (14), dimensionless", g.gammaIsasr);""","""    cmd.AddValue("gamma", "isasr: gamma of eq. (14), dimensionless", g.gammaIsasr);
    cmd.AddValue("isasrIslAct", "isasr: activeness cost on inter-satellite links too, as in Algorithm 3", g.isasrIslAct);""")
open('ntn-constellation-route-lifetime.cc','w').write(s)
# ---- router: cost-based Dijkstra
h=open('contact-graph-router.h').read()
a='''    /**
     * \\brief Longest-lived path within a delay stretch of the shortest one.'''
assert h.count(a)==1
h=h.replace(a,'''    /// P7: Dijkstra with an edge cost computed from (u, v, weight); returned total_weight is the true weight sum.
    using EdgeCost = std::function<double(uint32_t, uint32_t, double)>;
    WeightedPath ShortestPathCost(uint32_t src, uint32_t dst, const EdgeCost& cost) const;

'''+a)
open('contact-graph-router.h','w').write(h)
c=open('contact-graph-router.cc').read()
a='std::vector<uint32_t>\nContactGraphRouter::ShortestPath(uint32_t src, uint32_t dst) const'
assert c.count(a)==1
c=c.replace(a,'''ContactGraphRouter::WeightedPath
ContactGraphRouter::ShortestPathCost(uint32_t src, uint32_t dst, const EdgeCost& cost) const
{
    ++m_queries;
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
    using QEntry = std::pair<double, uint32_t>;
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
            double total = 0.0;
            while (c != src)
            {
                total += EdgeWeight(c, m_pred[c]);
                c = m_pred[c];
                path.push_back(c);
            }
            std::reverse(path.begin(), path.end());
            return {std::move(path), total};
        }
        if (w > m_dist[u])
        {
            continue;
        }
        for (const auto& [v, we] : m_flat[u])
        {
            const double ce = cost(u, v, we);
            if (!std::isfinite(ce))
            {
                continue;
            }
            const double newW = w + ce;
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

'''+a)
open('contact-graph-router.cc','w').write(c)
