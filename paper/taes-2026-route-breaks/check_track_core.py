#!/usr/bin/env python3
"""Exact-model Monte Carlo of the short-arc residual estimate (Lemma 2): circular Walker orbits under a rotating
spherical Earth with the gateway on the WGS84 ellipsoid, as in the simulator up to J2. Random visible satellite and time;
the true residual to the drop elevation against the parabola through psi^2 at t, t-15 s, t-30 s. Writes out/track_check.json."""
import json, numpy as np, isl_core as C
SH = {"starlink": (72, 22, 39, 53.0, 550e3, False), "kuiper": (34, 34, 1, 51.9, 630e3, False), "oneweb": (12, 49, 1, 87.9, 1200e3, True)}
GW = {"London": (51.51, -0.13), "New York": (40.71, -74.01), "Sydney": (-33.87, 151.21), "Los Angeles": (34.05, -118.24)}
WE = 7.2921159e-5; A, E2 = 6378137.0, 6.69437999014e-3; DROP = 23.0; UP = 25.0


def gw(lat, lon):
    la, lo = np.radians(lat), np.radians(lon); nn = A / np.sqrt(1 - E2 * np.sin(la) ** 2)
    return np.array([nn * np.cos(la) * np.cos(lo), nn * np.cos(la) * np.sin(lo), nn * (1 - E2) * np.sin(la)]), \
        np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def elev(sh, k, t, g, up):
    w = np.sqrt(C.MU / sh["r"] ** 3)
    u = sh["u0"][k] + w * t; O = sh["Om"][k] - WE * t; i = sh["i"]; r = sh["r"]
    x = r * np.stack([np.cos(O) * np.cos(u) - np.sin(O) * np.sin(u) * np.cos(i), np.sin(O) * np.cos(u) + np.cos(O) * np.sin(u) * np.cos(i), np.sin(u) * np.sin(i)], -1)
    d = x - g
    return np.degrees(np.arcsin((d @ up) / np.linalg.norm(d, axis=-1)))


def psi(e, r):
    e = np.radians(e); return np.arccos(6371000.0 * np.cos(e) / r) - e


