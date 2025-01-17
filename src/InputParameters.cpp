#include "InputParameters.h"

#include <algorithm>
#include <iostream>
#include <cmath>

#include "pugixml.hpp"

void
parseVecFromString(const std::string & str, std::vector<double> & res)
{
  auto s = str;
  s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
  s.erase(std::remove(s.begin(), s.end(), '\n'), s.end());

  auto current_delim_pos = s.find(",");
  std::size_t offset = 0u;
  auto previous_delim_pos = 0u;
  while (current_delim_pos != std::string::npos)
  {
    res.emplace_back(std::stod(s.substr(
        previous_delim_pos + offset, current_delim_pos - (previous_delim_pos + offset))));
    previous_delim_pos = current_delim_pos;
    current_delim_pos = s.find(",", current_delim_pos + 1);
    offset = 1u;
  }
  res.emplace_back(std::stod(s.substr(previous_delim_pos + offset)));
}

void
parseVecFromString(const std::string & str, std::vector<unsigned int> & res)
{
  auto s = str;
  s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
  s.erase(std::remove(s.begin(), s.end(), '\n'), s.end());

  auto current_delim_pos = s.find(",");
  std::size_t offset = 0u;
  auto previous_delim_pos = 0u;
  while (current_delim_pos != std::string::npos)
  {
    res.emplace_back(std::stoul(s.substr(
        previous_delim_pos + offset, current_delim_pos - (previous_delim_pos + offset))));
    previous_delim_pos = current_delim_pos;
    current_delim_pos = s.find(",", current_delim_pos + 1);
    offset = 1u;
  }
  res.emplace_back(std::stoul(s.substr(previous_delim_pos + offset)));
}

double
parseDoubleParam(const pugi::xml_node & node, const std::string & p_name)
{
  if (node.attribute(p_name.c_str()).as_double(-1.0) >= 0.0)
    return node.attribute(p_name.c_str()).as_double(-1.0);
  else
  {
    std::cerr << "You must provide a valid param in " << p_name << "!" << std::endl;
    std::exit(1);
    return -1.0;
  }
}

unsigned int
parseUintParam(const pugi::xml_node & node, const std::string & p_name)
{
  if (node.attribute(p_name.c_str()).as_uint(std::numeric_limits<unsigned int>::max()) != std::numeric_limits<unsigned int>::max())
    return node.attribute(p_name.c_str()).as_uint(0);
  else
  {
    std::cerr << "You must provide a valid param in " << p_name << "!" << std::endl;
    std::exit(1);
    return 0;
  }
}

std::string
parseStringParam(const pugi::xml_node & node, const std::string & p_name)
{
  if (std::string(node.attribute(p_name.c_str()).as_string("")) != "")
    return std::string(node.attribute(p_name.c_str()).as_string(""));
  else
  {
    std::cerr << "You must provide a valid param in " << p_name << "!" << std::endl;
    std::exit(1);
    return "";
  }
}

// I hate input file parsing...
InputParameters
parseInputParameters(const std::string & file_path)
{
  pugi::xml_document doc;
  auto result = doc.load_file(file_path.c_str());
  if (!result)
  {
    std::cerr << "Failed to parse " << file_path << std::endl;
    std::cerr << result.description() << std::endl;
    std::exit(1);
  }

  auto p = doc.child("doma_inp");
  InputParameters params;

  // Parse solver input parameters.
  {
    auto s = p.child("solver");
    if (!s)
    {
      std::cerr << "A solver block has not been provided!" << std::endl;
      std::exit(1);
    }

    if (std::string(s.attribute("solve_type").as_string()) == "steady")
      params._mode = RunMode::FixedSrc;
    else if(std::string(s.attribute("solve_type").as_string()) == "eigen")
      params._mode = RunMode::Eigen;
    else if(std::string(s.attribute("solve_type").as_string()) == "transient")
      params._mode = RunMode::Transient;
    else
    {
      std::cerr << "Unknown run_mode " << s.attribute("solve_type").as_string() << std::endl;
      std::exit(1);
    }

    if (std::string(s.attribute("scheme").as_string()) == "dd")
      params._eq_type = EquationType::DD;
    else if(std::string(s.attribute("scheme").as_string()) == "tw_dd")
      params._eq_type = EquationType::TW_DD;
    else
    {
      std::cerr << "Unknown scheme " << s.attribute("scheme").as_string() << std::endl;
      std::exit(1);
    }

    params._num_e_groups = parseUintParam(s, "num_groups");

    params._src_it_tol = parseDoubleParam(s, "fs_tol");
    params._num_src_it = parseUintParam(s, "fs_max_it");
    params._gs_tol = parseDoubleParam(s, "mg_tol");
    params._num_mg_it = parseUintParam(s, "mg_max_it");

    if (params._mode == RunMode::Eigen)
      params._k_tol = parseDoubleParam(s, "k_tol");

    if (params._mode == RunMode::Transient)
    {
      params._t0 = parseDoubleParam(s, "t0");
      params._t1 = parseDoubleParam(s, "t1");
      params._num_steps = parseUintParam(s, "t_steps");

      if (std::string(s.attribute("ics").as_string()) == "zero")
        params._ic = TransientIC::Zero;
      else if(std::string(s.attribute("ics").as_string()) == "steady")
        params._ic = TransientIC::SteadyState;
      else
      {
        std::cerr << "Unknown initial condition " << s.attribute("ics").as_string() << std::endl;
        std::exit(1);
      }
    }
  }
  // End parsing solver parameters.

  // Parsing the angular quadrature parameters.
  {
    auto aq = p.child("angular_quad");
    if (!aq)
    {
      std::cerr << "An angular quadrature block has not been provided!" << std::endl;
      std::exit(1);
    }
    params._num_polar = parseUintParam(aq, "np");
    params._num_azimuthal = parseUintParam(aq, "na");
    // type doesn't matter for now. TODO: add other angular quadrature rules.
  }
  // End parsing the angular quadrature parameters

  // Parsing the mesh.
  {
    auto m = p.child("mesh");
    if (!m)
    {
      std::cerr << "A mesh block has not been provided!" << std::endl;
      std::exit(1);
    }

    params._num_dims = parseUintParam(m, "dims");
    params._refinement = parseUintParam(m, "refinement");
    {
      auto dx = m.child("dx");
      if (!dx)
      {
        std::cerr << "Mesh x dimensions have not been provided!" << std::endl;
        std::exit(1);
      }
      parseVecFromString(parseStringParam(dx, "vals"), params._dx);

      auto ix = m.child("ix");
      if (!ix)
      {
        std::cerr << "Mesh x divisions have not been provided!" << std::endl;
        std::exit(1);
      }
      parseVecFromString(parseStringParam(ix, "vals"), params._x_intervals);
      for (auto & i : params._x_intervals)
        i *= (params._refinement + 1);
    }

    if (params._num_dims > 1)
    {
      auto dy = m.child("dy");
      if (!dy)
      {
        std::cerr << "Mesh y dimensions have not been provided!" << std::endl;
        std::exit(1);
      }
      parseVecFromString(parseStringParam(dy, "vals"), params._dy);

      auto iy = m.child("iy");
      if (!iy)
      {
        std::cerr << "Mesh y divisions have not been provided!" << std::endl;
        std::exit(1);
      }
      parseVecFromString(parseStringParam(iy, "vals"), params._y_intervals);
      for (auto & i : params._y_intervals)
        i *= (params._refinement + 1);
    }

    if (params._num_dims > 2)
    {
      auto dz = m.child("dz");
      if (!dz)
      {
        std::cerr << "Mesh z dimensions have not been provided!" << std::endl;
        std::exit(1);
      }
      parseVecFromString(parseStringParam(dz, "vals"), params._dz);

      auto iz = m.child("iz");
      if (!iz)
      {
        std::cerr << "Mesh z divisions have not been provided!" << std::endl;
        std::exit(1);
      }
      parseVecFromString(parseStringParam(iz, "vals"), params._z_intervals);
      for (auto & i : params._z_intervals)
        i *= (params._refinement + 1);
    }

    {
      auto b = m.child("blocks");
      if (!b)
      {
        std::cerr << "Mesh blocks have not been provided!" << std::endl;
        std::exit(1);
      }
      parseVecFromString(parseStringParam(b, "vals"), params._blocks);
    }
  }
  // End parsing the mesh parameters.

  // Parsing material properties.
  {
    auto mat_block = p.child("materials");
    if (!mat_block)
    {
      std::cerr << "A materials block has not been provided!" << std::endl;
      std::exit(1);
    }

    for (auto & mat_node : mat_block)
    {
      auto blk = parseUintParam(mat_node, "block");
      if (params._block_mat_info.count(blk) != 0u)
      {
        std::cerr << "A block can only contain one set of material properties, "
                  << "you have provided multiple material properties for block "
                  << blk << "!" << std::endl;
        std::exit(1);
      }

      auto & mat_props = (*params._block_mat_info.emplace(blk, MaterialProps()).first).second;
      if (params._mode == RunMode::Transient)
        mat_props._num_d_groups = parseUintParam(mat_node, "delayed_groups");
      else
        mat_props._num_d_groups = 0u;

      for (auto & rxn : mat_node)
      {
        if (std::string(rxn.attribute("type").as_string()) == "total")
          parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._g_total);
        else if (std::string(rxn.attribute("type").as_string()) == "scatter")
        {
          parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._g_g_scatter_mat);
          if (mat_props._g_g_scatter_mat.size() != (params._num_e_groups * params._num_e_groups))
          {
            std::cerr << "Not enough cross sections have been provided! The simulation requires "
                      << (params._num_e_groups * params._num_e_groups) << " cross sections; "
                      << "'scatter' only provides " << mat_props._g_g_scatter_mat.size() << "." << std::endl;
            std::exit(1);
          }

          if (rxn.attribute("transpose").as_bool(false))
          {
            auto tr = mat_props._g_g_scatter_mat;
            for (unsigned int r = 0u; r < params._num_e_groups; ++r)
              for (unsigned int c = 0u; c < params._num_e_groups; ++c)
                mat_props._g_g_scatter_mat[r * params._num_e_groups + c] = tr[c * params._num_e_groups + r];
          }
        }
        else if (std::string(rxn.attribute("type").as_string()) == "source")
           parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._g_src);
        else if (std::string(rxn.attribute("type").as_string()) == "chi_p")
           parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._g_chi_p);
        else if (std::string(rxn.attribute("type").as_string()) == "nu_sigma_f")
           parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._g_prod);
        else if (std::string(rxn.attribute("type").as_string()) == "inv_vel")
          parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._g_inv_v);
        else if (std::string(rxn.attribute("type").as_string()) == "chi_d")
          parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._n_g_chi_d);
        else if (std::string(rxn.attribute("type").as_string()) == "beta")
          parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._g_n_beta);
        else if (std::string(rxn.attribute("type").as_string()) == "lambda_d")
          parseVecFromString(parseStringParam(rxn, "mgxs"), mat_props._n_lambda);
        else
        {
          std::cerr << "Unsupported reaction type " << rxn.attribute("type").as_string() << std::endl;
          std::exit(1);
        }
      }

      if (mat_props._g_total.size() != params._num_e_groups)
      {
        std::cerr << "Not enough cross sections have been provided! The simulation requires "
                  << params._num_e_groups << " cross sections; 'total' only provides "
                  << mat_props._g_total.size() << "." << std::endl;
        std::exit(1);
      }
      if (mat_props._g_g_scatter_mat.size() != (params._num_e_groups * params._num_e_groups))
      {
        std::cerr << "Not enough cross sections have been provided! The simulation requires "
                  << (params._num_e_groups * params._num_e_groups) << " cross sections; "
                  << "'scatter' only provides " << mat_props._g_g_scatter_mat.size() << "." << std::endl;
        std::exit(1);
      }
      if (mat_props._g_src.size() != params._num_e_groups && mat_props._g_src.size() != 0u)
      {
        std::cerr << "Not enough external sources have been provided! The simulation requires "
                  << params._num_e_groups << " sources; 'source' only provides "
                  << mat_props._g_src.size() << "." << std::endl;
        std::exit(1);
      }
      if (mat_props._g_chi_p.size() != params._num_e_groups && mat_props._g_chi_p.size() != 0u)
      {
        std::cerr << "Not enough prompt fission spectra values have been provided! The simulation requires "
                  << params._num_e_groups << " prompt fission spectra values; 'chi' only provides "
                  << mat_props._g_chi_p.size() << "." << std::endl;
        std::exit(1);
      }
      if (mat_props._g_prod.size() != params._num_e_groups && mat_props._g_prod.size() != 0u)
      {
        std::cerr << "Not enough fission production cross sections have been provided! The simulation requires "
                  << params._num_e_groups << " fission production cross sections; 'nu_sigma_f' only provides "
                  << mat_props._g_prod.size() << "." << std::endl;
        std::exit(1);
      }
      if (mat_props._g_inv_v.size() != params._num_e_groups && params._mode == RunMode::Transient)
      {
        std::cerr << "Not enough inverse velocities have been provided! The simulation requires "
                  << params._num_e_groups << " inverse velocities when running in transient mode; "
                  << "'inv_vel' only provides " << mat_props._g_inv_v.size() << "." << std::endl;
        std::exit(1);
      }
      if (mat_props._n_g_chi_d.size() != params._num_e_groups * mat_props._num_d_groups && params._mode == RunMode::Transient)
      {
        std::cerr << "Not enough delayed fission spectra values have been provided! The simulation requires "
                  << params._num_e_groups * mat_props._num_d_groups << " delayed fission spectra values when "
                  "running in transient mode; 'chi_d' only provides " << mat_props._n_g_chi_d.size() << "." << std::endl;
        std::exit(1);
      }
      if (mat_props._g_n_beta.size() != params._num_e_groups * mat_props._num_d_groups && params._mode == RunMode::Transient)
      {
        std::cerr << "Not enough delayed fractions have been provided! The simulation requires "
                  << params._num_e_groups * mat_props._num_d_groups << " delayed fractions when "
                  "running in transient mode; 'beta' only provides " << mat_props._g_n_beta.size() << "." << std::endl;
        std::exit(1);
      }
      if (mat_props._n_lambda.size() != mat_props._num_d_groups && params._mode == RunMode::Transient)
      {
        std::cerr << "Not enough delayed decay constants have been provided! The simulation requires "
                  << mat_props._num_d_groups << " delayed decay constants when "
                  "running in transient mode; 'lambda_d' only provides " << mat_props._n_lambda.size() << "." << std::endl;
        std::exit(1);
      }
    }
  }
  // End parsing the material properties.

  // Parsing transients.
  {
    auto tr_block = p.child("transients");
    if (tr_block && params._mode == RunMode::Transient)
    {
      for (auto & tr_node : tr_block)
      {
        // This is a source step insertion / step removal transient.
        if (std::string(tr_node.name()) == "source_step")
        {
          auto blk = parseUintParam(tr_node, "block");
          if (params._block_step_src.count(blk) != 0u)
          {
            std::cerr << "A block can only contain one step source, "
                      << "you have provided multiple step sources for block "
                      << blk << "!" << std::endl;
            std::exit(1);
          }

          auto & step = (*params._block_step_src.emplace(blk, SourceStep()).first).second;
          if (std::string(tr_node.attribute("insert_at").as_string()) != "" &&
              std::string(tr_node.attribute("remove_at").as_string()) != "")
          {
            step._type = StepType::Both;
            step._insert_time = parseDoubleParam(tr_node, "insert_at");
            step._remove_time = parseDoubleParam(tr_node, "remove_at");
            if (step._remove_time < step._insert_time)
            {
              std::cerr << "The step source for block " << blk
                        << " is being removed at " << step._remove_time
                        << ", which is less than the insertion time of "
                        << step._insert_time << "!" << std::endl;
              std::exit(1);
            }

            if (step._insert_time < params._t0 || step._remove_time < params._t0)
            {
              std::cerr << "The step source for block " << blk
                        << " is being inserted / removed before "
                        << "the simulation starts!" << std::endl;
            }
            if (step._insert_time > params._t1 || step._remove_time > params._t1)
            {
              std::cerr << "The step source for block " << blk
                        << " is being inserted / removed after "
                        << "the simulation ends!" << std::endl;
            }
          }
          else if (std::string(tr_node.attribute("insert_at").as_string()) != "")
          {
            step._type = StepType::Insert;
            step._insert_time = parseDoubleParam(tr_node, "insert_at");

            if (step._insert_time < params._t0)
            {
              std::cerr << "The step source for block " << blk
                        << " is being inserted before "
                        << "the simulation starts!" << std::endl;
            }
            if (step._insert_time > params._t1)
            {
              std::cerr << "The step source for block " << blk
                        << " is being inserted after "
                        << "the simulation ends!" << std::endl;
            }
          }
          else if (std::string(tr_node.attribute("remove_at").as_string()) != "")
          {
            step._type = StepType::Remove;
            step._remove_time = parseDoubleParam(tr_node, "remove_at");

            if (step._remove_time < params._t0)
            {
              std::cerr << "The step source for block " << blk
                        << " is being removed before "
                        << "the simulation starts!" << std::endl;
            }
            if (step._remove_time > params._t1)
            {
              std::cerr << "The step source for block " << blk
                        << " is being removed after "
                        << "the simulation ends!" << std::endl;
            }
          }
          else
          {
            std::cerr << "Invalid step source for block " << blk << "!" << std::endl;
            std::exit(1);
          }

          parseVecFromString(parseStringParam(tr_node, "mgxs"), step._g_src);
          if (step._g_src.size() != params._num_e_groups)
          {
            std::cerr << "Not enough external sources have been provided for the step source in block " << blk
                      << "! The simulation requires " << params._num_e_groups << " sources; "
                      << step._g_src.size() << " have been provided." << std::endl;
            std::exit(1);
          }
        }
      }
    }
  }
  // End parsing transients.

  return params;
}
