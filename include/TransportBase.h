#pragma once

enum class CertesianAxis
{
  X = 0u,
  Y = 1u,
  Z = 2u
};

enum class CertesianFaceSide
{
  Front = 0u, // Face normal along +y.
  Back = 1u,  // Face normal along -y.
  Right = 2u, // Face normal along +x.
  Left = 3u,  // Face normal along -x.
  Top = 4u,   // Face normal along +z.
  Bottom = 5u // Face normal along -z.
};

enum class Octant
{
  PPP = 0u, // +\mu, +\eta, +\xi
  PPM = 1u, // +\mu, +\eta, -\xi
  PMP = 2u, // +\mu, -\eta, +\xi
  PMM = 3u, // +\mu, -\eta, -\xi
  MPP = 4u, // -\mu, +\eta, +\xi
  MPM = 5u, // -\mu, +\eta, -\xi
  MMP = 6u, // -\mu, -\eta, +\xi
  MMM = 7u, // -\mu, -\eta, -\xi
};

enum class DiscretizationType
{
  DiamondDifference = 0u,
  StepCharacteristics = 1u
};
