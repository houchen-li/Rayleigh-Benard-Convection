#include <utility>

#include <utility>

#include <iostream>
#include <boost/program_options.hpp>
#include <H5Cpp.h>
#include "../include/def.h"
#include "../include/utils.h"
#include "../include/state_type.h"
#include "../include/observer_type.h"
#include "../include/solver_type.h"
#include "../include/herring_func.hpp"


class HerringSolverType: public RBCSystem::SolverType {
public:
    explicit HerringSolverType(const Real &Prandtl_number = 1.0, const Real &Rayleigh_number = 1000.0,
                               const Real &wave_number = 1.0, const Real &start_time = 0.0, const Real &end_time = 1.0,
                               const String &in_h5file_name = static_cast<String>("in_data.h5"),
                               const String &in_h5dset_name = static_cast<String>("/group/state"),
                               const String &out_h5file_name = static_cast<String>("out_data.h5"),
                               const String &out_h5group_name = static_cast<String>("/group/state"));

    HerringSolverType(const HerringSolverType &rhs);

    HerringSolverType &operator=(const HerringSolverType &rhs);

    ~HerringSolverType(void);

    const String &in_h5file_name(void) const;

    void in_h5file_name(const String &in_h5file_name);

    const String &in_h5dset_name(void) const;

    void in_h5dset_name(const String &in_h5dset_name);

    const String &out_h5file_name(void) const;

    void out_h5file_name(const String &out_h5file_name);

    const String &out_h5group_name(void) const;

    void out_h5group_name(const String &out_h5group_name);

    void operator()(const RBCSystem::StateType &f, RBCSystem::StateType &dfdt, const Real &t) const override;

    void parseArgv(int argc, char *argv[]);

private:
    String inf_name, ind_name, outf_name, outg_name;
};


int main(int argc, char* argv[])
{
    HerringSolverType solver;
    H5::H5File inf, outf;

    solver.parseArgv(argc, argv);
    //inf = H5::H5File(mysolver.in_h5file_name, H5F_ACC_RDONLY);
    //initial_state.loadFile(inf, mysolver.in_h5obj_name);
    //inf.close();
    //solver.initial_state(initial_state);
    solver.evaluate();
    outf = H5::H5File(solver.out_h5file_name(), H5F_ACC_TRUNC);
    solver.saveFile(outf, solver.out_h5group_name());
    outf.close();

    return 0;
}

HerringSolverType::HerringSolverType(const Real &Prandtl_number, const Real &Rayleigh_number, const Real &wave_number,
        const Real &start_time, const Real &end_time,
        const String &in_h5file_name, const String &in_h5dset_name,
        const String &out_h5file_name, const String &out_h5group_name):
        RBCSystem::SolverType(Prandtl_number, Rayleigh_number, wave_number, RBCSystem::StateType::TrivialState(), start_time, end_time),
        inf_name(in_h5file_name), ind_name(in_h5dset_name),
        outf_name(out_h5file_name), outg_name(out_h5group_name) { return; }

HerringSolverType::HerringSolverType(const HerringSolverType &rhs):
        RBCSystem::SolverType(rhs),
        inf_name(rhs.inf_name), ind_name(rhs.ind_name),
        outf_name(rhs.inf_name), outg_name(rhs.outg_name) { return; }

HerringSolverType &HerringSolverType::operator=(const HerringSolverType &rhs)
{
    RBCSystem::SolverType::operator=(rhs);
    inf_name = rhs.inf_name;
    ind_name = rhs.ind_name;
    outf_name = rhs.outf_name;
    outg_name = rhs.outg_name;

    return *this;
}

HerringSolverType::~HerringSolverType(void) { return; }

inline const String &HerringSolverType::in_h5file_name(void) const { return inf_name; }

inline void HerringSolverType::in_h5file_name(const String &in_h5file_name) { inf_name = in_h5file_name; return; }

inline const String &HerringSolverType::in_h5dset_name(void) const { return ind_name; }

inline void HerringSolverType::in_h5dset_name(const String &in_h5dset_name) { ind_name = in_h5dset_name; return; }

inline const String &HerringSolverType::out_h5file_name(void) const { return outf_name; }

inline void HerringSolverType::out_h5file_name(const String &out_h5file_name) { outf_name = out_h5file_name; return; }

inline const String &HerringSolverType::out_h5group_name(void) const { return outg_name; }

inline void HerringSolverType::out_h5group_name(const String &out_h5group_name) { outg_name = out_h5group_name; return; }

void HerringSolverType::operator()(const RBCSystem::StateType &f, RBCSystem::StateType &dfdt, const Real &t) const
{
    func(f, dfdt, t, Prandtl_number(), Rayleigh_number(), wave_number());

    return;
}

void HerringSolverType::parseArgv(int argc, char *argv[])
{
    using namespace boost::program_options;

    options_description desc("**This program is used to simulate HerringSolver's model.**\nUsages:");
    desc.add_options()
            ("help,h", "display usages;")
            ("Prandtl_number", value<Real>()->default_value(1.0), "Prandl_number;")
            ("Rayleigh_number", value<Real>()->default_value(1000.0), "Rayleigh_number;")
            ("wave_number", value<Real>()->default_value(0.5), "wave_number;")
            ("start_time", value<Real>()->default_value(0.0), "start_time")
            ("end_time", value<Real>()->default_value(1.0), "end_time")
            ("input_h5file", value<String>()->default_value(static_cast<String>("in_data.h5")), "input h5file's name;")
            ("input_h5obj", value<String>()->default_value(static_cast<String>("/group/state")), "input h5obj's name;")
            ("output_h5file", value<String>()->default_value(static_cast<String>("out_data.h5")), "output h5file's name;")
            ("output_h5obj", value<String>()->default_value(static_cast<String>("/group")), "output h5obj's name.");
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if(vm.count("help") ){
        std::cout << desc << std::endl;
        std::exit(1);
    }
    Prandtl_number(vm["Prandtl_number"].as<Real>());
    Rayleigh_number(vm["Rayleigh_number"].as<Real>());
    wave_number(vm["wave_number"].as<Real>());
    start_time(vm["start_time"].as<Real>());
    end_time(vm["end_time"].as<Real>());
    in_h5file_name(vm["input_h5file"].as<String>());
    in_h5dset_name(vm["input_h5obj"].as<String>());
    out_h5file_name(vm["output_h5file"].as<String>());
    out_h5group_name(vm["output_h5obj"].as<String>());

    return;
}
