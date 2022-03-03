/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "atoms.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <type_traits>
#include <vector>

#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

const char* enumValueToString(ParticleType enumValue)
{
    static constexpr gmx::EnumerationArray<ParticleType, const char*> particleTypeNames = {
        "Atom", "Nucleus", "Shell", "Bond", "VSite"
    };
    return particleTypeNames[enumValue];
}

namespace
{

/*! \brief Serialize unsigned short values.
 *
 * \tparam T type to be serialized
 * \param[in,out] serializer the serializer
 * \param[in,out] value to be serialized
 */
template<typename T>
std::enable_if_t<std::is_same_v<T, unsigned short>, void> serializeValue(gmx::ISerializer* serializer,
                                                                         T*                value)
{
    serializer->doUShort(value);
}

/*! \brief Serialize real values.
 *
 * \tparam T type to be serialized
 * \param[in,out] serializer the serializer
 * \param[in,out] value to be serialized
 */
template<typename T>
std::enable_if_t<std::is_same_v<T, real>, void> serializeValue(gmx::ISerializer* serializer, T* value)
{
    serializer->doReal(value);
}

/*! \brief Serialize NameHolder values.
 *
 * \tparam T type to be serialized
 * \param[in,out] serializer the serializer
 * \param[in,out] value to be serialized
 */
template<typename T>
std::enable_if_t<std::is_same_v<T, NameHolder>, void> serializeValue(gmx::ISerializer* serializer, T* value)
{
    GMX_ASSERT(!serializer->reading(), "This only works with a writing serializer");
    GMX_ASSERT(value->has_value(), "Can not access uninitialized element");
    (*value)->serialize(serializer);
}

/*! \brief Deserialize NameHolder values.
 *
 * \tparam T type to be deerialized
 * \param[in,out] serializer the serializer
 * \param[in,out] value to be serialized
 */
template<typename T>
std::enable_if_t<std::is_same_v<T, NameHolder>, void> serializeValue(gmx::ISerializer*  serializer,
                                                                     T*                 value,
                                                                     const StringTable& table)
{
    GMX_ASSERT(serializer->reading(), "This only works with a reading serializer");
    *value = readStringTableEntry(serializer, table);
}
} // namespace

template<typename T>
void FEPStateValue<T>::serialize(gmx::ISerializer* serializer)
{
    GMX_ASSERT(!serializer->reading(), "Can not write with reading serializer");
    serializer->doBool(&haveBState_);
    serializeValue<T>(serializer, &storage_[0]);
    if (haveBState_)
    {
        serializeValue<T>(serializer, &storage_[1]);
    }
}

template<typename T>
FEPStateValue<T>::FEPStateValue(gmx::ISerializer* serializer)
{
    GMX_ASSERT(serializer->reading(), "Can not create with writing serializer");
    serializer->doBool(&haveBState_);
    serializeValue<T>(serializer, &storage_[0]);
    if (haveBState_)
    {
        serializeValue<T>(serializer, &storage_[1]);
    }
}

template<typename T>
FEPStateValue<T>::FEPStateValue(gmx::ISerializer* serializer, const StringTable& table)
{
    GMX_ASSERT(serializer->reading(), "Can not create with writing serializer");
    serializer->doBool(&haveBState_);
    serializeValue<T>(serializer, &storage_[0], table);
    if (haveBState_)
    {
        serializeValue<T>(serializer, &storage_[1], table);
    }
}

SimulationParticle::SimulationParticle(gmx::ISerializer* serializer, const StringTable& table)
{
    GMX_ASSERT(serializer->reading(), "Can not create particle with writing serializer");
    mass_              = ParticleMass(serializer);
    charge_            = ParticleCharge(serializer);
    particleTypeValue_ = ParticleTypeValue(serializer);
    particleTypeName_  = ParticleTypeName(serializer, table);
    particleName_      = readStringTableEntry(serializer, table);
    serializer->doEnumAsInt<ParticleType>(&particleType_);
    int residueIndex;
    serializer->doInt(&residueIndex);
    residueIndex_ = residueIndex;
    serializer->doInt(&atomicNumber_);
    serializer->doBool(&haveMass_);
    serializer->doBool(&haveCharge_);
    serializer->doBool(&haveType_);
    serializer->doBool(&haveParticleName_);
    serializer->doBool(&haveParticleTypeName_);
    haveBState_ = mass_.haveBState_ && charge_.haveBState_ && particleTypeValue_.haveBState_
                  && particleTypeName_.haveBState_;
}

void SimulationParticle::serializeParticle(gmx::ISerializer* serializer)
{
    GMX_ASSERT(!serializer->reading(), "Can not write particle with reading serializer");
    mass_.serialize(serializer);
    charge_.serialize(serializer);
    particleTypeValue_.serialize(serializer);
    particleTypeName_.serialize(serializer);
    particleName_->serialize(serializer);
    serializer->doEnumAsInt<ParticleType>(&particleType_);
    int residueIndex = residueIndex_;
    serializer->doInt(&residueIndex);
    serializer->doInt(&atomicNumber_);
    serializer->doBool(&haveMass_);
    serializer->doBool(&haveCharge_);
    serializer->doBool(&haveType_);
    serializer->doBool(&haveParticleName_);
    serializer->doBool(&haveParticleTypeName_);
}

SimulationResidue::SimulationResidue(gmx::ISerializer* serializer, const StringTable& table)
{
    GMX_ASSERT(serializer->reading(), "Can not create residue with writing serializer");
    residueName_ = readStringTableEntry(serializer, table);
    serializer->doInt64(&nr_);
    serializer->doUChar(&insertionCode_);
    serializer->doInt64(&begin_);
    serializer->doInt64(&size_);
}

void SimulationResidue::serializeResidue(gmx::ISerializer* serializer)
{
    GMX_ASSERT(!serializer->reading(), "Can not write residue with reading serializer");
    residueName_->serialize(serializer);
    serializer->doInt64(&nr_);
    serializer->doUChar(&insertionCode_);
    serializer->doInt64(&begin_);
    serializer->doInt64(&size_);
}

real PdbEntry::occup() const
{
    if (!hasOccupancy_)
    {
        GMX_THROW(gmx::InternalError("Can not access occupancy information that is not present"));
    }
    else
    {
        return occupancy_;
    }
}

real PdbEntry::bfac() const
{
    if (!hasbFactor_)
    {
        GMX_THROW(gmx::InternalError("Can not access B-Factor information that is not present"));
    }
    else
    {
        return bFactor_;
    }
}

gmx::ArrayRef<const int> PdbEntry::uij() const
{
    if (!hasAnistropy_)
    {
        GMX_THROW(gmx::InternalError("Can not access anisotropy information that is not present"));
    }
    else
    {
        return uij_;
    }
}

SimulationResidue SimulationResidueBuilder::finalize(gmx::index residueNumber,
                                                     gmx::index begin,
                                                     gmx::index end) const
{
    return SimulationResidue(name_, residueNumber, insertionCode_, begin, end);
}

SimulationMolecule SimulationMoleculeBuilder::finalize()
{
    std::vector<SimulationResidue> residues;
    int                            residueNumber = 0;
    for (const auto& builderResidue : residues_)
    {
        const int first = particles_.size();
        for (const auto& particle : builderResidue.particles())
        {
            addParticle(particle, residueNumber);
        }
        const int size = particles_.size() - first;

        residues.emplace_back(builderResidue.finalize(residueNumber, first, size));
        residueNumber++;
    }

    SimulationMolecule molecule(
            &particles_,
            &residues,
            std::all_of(particles_.begin(),
                        particles_.end(),
                        [](const auto& particle) { return particle.haveMass(); }),
            std::all_of(particles_.begin(),
                        particles_.end(),
                        [](const auto& particle) { return particle.haveCharge(); }),
            std::all_of(particles_.begin(),
                        particles_.end(),
                        [](const auto& particle) { return particle.haveParticleName(); }),
            std::all_of(particles_.begin(),
                        particles_.end(),
                        [](const auto& particle) { return particle.haveType(); }),
            std::all_of(particles_.begin(),
                        particles_.end(),
                        [](const auto& particle) { return particle.haveParticleTypeName(); }),
            std::all_of(particles_.begin(), particles_.end(), [](const auto& particle) {
                return particle.haveBState();
            }));
    return molecule;
}

SimulationMolecule::SimulationMolecule(std::vector<SimulationParticle>* particles,
                                       std::vector<SimulationResidue>*  residues,
                                       bool                             allParticlesHaveMass,
                                       bool                             allParticlesHaveCharge,
                                       bool                             allParticlesHaveAtomName,
                                       bool                             allParticlesHaveType,
                                       bool                             allParticlesHaveTypeName,
                                       bool                             allParticlesHaveBstate) :
    allParticlesHaveMass_(allParticlesHaveMass),
    allParticlesHaveCharge_(allParticlesHaveCharge),
    allParticlesHaveAtomName_(allParticlesHaveAtomName),
    allParticlesHaveType_(allParticlesHaveType),
    allParticlesHaveTypeName_(allParticlesHaveTypeName),
    allParticlesHaveBstate_(allParticlesHaveBstate)
{
    std::swap(*particles, particles_);
    std::swap(*residues, residues_);
}

SimulationMolecule::SimulationMolecule(gmx::ISerializer* serializer, const StringTable& table)
{
    GMX_ASSERT(serializer->reading(), "Can not read molecule with writing serializer");
    int numParticle;
    serializer->doInt(&numParticle);
    for (int index = 0; index < numParticle; ++index)
    {
        particles_.emplace_back(SimulationParticle(serializer, table));
    }
    int numResidue;
    serializer->doInt(&numResidue);
    for (int index = 0; index < numResidue; ++index)
    {
        residues_.emplace_back(SimulationResidue(serializer, table));
    }
    serializer->doBool(&allParticlesHaveMass_);
    serializer->doBool(&allParticlesHaveCharge_);
    serializer->doBool(&allParticlesHaveAtomName_);
    serializer->doBool(&allParticlesHaveType_);
    serializer->doBool(&allParticlesHaveTypeName_);
    serializer->doBool(&allParticlesHaveBstate_);
}
void SimulationMolecule::serializeMolecule(gmx::ISerializer* serializer)
{
    GMX_ASSERT(!serializer->reading(), "Can not write molecule with reading serializer");
    int numParticle = particles_.size();
    serializer->doInt(&numParticle);
    for (auto& particle : particles_)
    {
        particle.serializeParticle(serializer);
    }
    int numResidue = residues_.size();
    serializer->doInt(&numResidue);
    for (auto& residue : residues_)
    {
        residue.serializeResidue(serializer);
    }
    serializer->doBool(&allParticlesHaveMass_);
    serializer->doBool(&allParticlesHaveCharge_);
    serializer->doBool(&allParticlesHaveAtomName_);
    serializer->doBool(&allParticlesHaveType_);
    serializer->doBool(&allParticlesHaveTypeName_);
    serializer->doBool(&allParticlesHaveBstate_);
}

void SimulationMoleculeBuilder::addParticle(const SimulationParticle& particle, gmx::index residueNumber)
{
    particles_.emplace_back(particle, residueNumber);
}

void SimulationMoleculeBuilder::addResidue(const SimulationResidueBuilder& residue)
{
    residues_.emplace_back(residue);
}

void init_atom(t_atoms* at)
{
    at->nr          = 0;
    at->nres        = 0;
    at->atom        = nullptr;
    at->resinfo     = nullptr;
    at->atomname    = nullptr;
    at->atomtype    = nullptr;
    at->atomtypeB   = nullptr;
    at->pdbinfo     = nullptr;
    at->haveMass    = FALSE;
    at->haveCharge  = FALSE;
    at->haveType    = FALSE;
    at->haveBState  = FALSE;
    at->havePdbInfo = FALSE;
}

void done_atom(t_atoms* at)
{
    sfree(at->atom);
    sfree(at->resinfo);
    sfree(at->atomname);
    sfree(at->atomtype);
    sfree(at->atomtypeB);
    sfree(at->pdbinfo);
    init_atom(at);
}

void done_and_delete_atoms(t_atoms* atoms)
{
    done_atom(atoms);
    delete atoms;
}

void add_t_atoms(t_atoms* atoms, int natom_extra, int nres_extra)
{
    if (natom_extra > 0)
    {
        srenew(atoms->atomname, atoms->nr + natom_extra);
        srenew(atoms->atom, atoms->nr + natom_extra);
        if (nullptr != atoms->pdbinfo)
        {
            srenew(atoms->pdbinfo, atoms->nr + natom_extra);
        }
        if (nullptr != atoms->atomtype)
        {
            srenew(atoms->atomtype, atoms->nr + natom_extra);
        }
        if (nullptr != atoms->atomtypeB)
        {
            srenew(atoms->atomtypeB, atoms->nr + natom_extra);
        }
        for (int i = atoms->nr; (i < atoms->nr + natom_extra); i++)
        {
            atoms->atomname[i] = nullptr;
            memset(&atoms->atom[i], 0, sizeof(atoms->atom[i]));
            if (nullptr != atoms->pdbinfo)
            {
                std::memset(&atoms->pdbinfo[i], 0, sizeof(atoms->pdbinfo[i]));
            }
            if (nullptr != atoms->atomtype)
            {
                atoms->atomtype[i] = nullptr;
            }
            if (nullptr != atoms->atomtypeB)
            {
                atoms->atomtypeB[i] = nullptr;
            }
        }
        atoms->nr += natom_extra;
    }
    if (nres_extra > 0)
    {
        srenew(atoms->resinfo, atoms->nres + nres_extra);
        for (int i = atoms->nres; (i < atoms->nres + nres_extra); i++)
        {
            std::memset(&atoms->resinfo[i], 0, sizeof(atoms->resinfo[i]));
        }
        atoms->nres += nres_extra;
    }
}

void init_t_atoms(t_atoms* atoms, int natoms, gmx_bool bPdbinfo)
{
    atoms->nr   = natoms;
    atoms->nres = 0;
    snew(atoms->atomname, natoms);
    atoms->atomtype  = nullptr;
    atoms->atomtypeB = nullptr;
    snew(atoms->resinfo, natoms);
    snew(atoms->atom, natoms);
    atoms->haveMass    = FALSE;
    atoms->haveCharge  = FALSE;
    atoms->haveType    = FALSE;
    atoms->haveBState  = FALSE;
    atoms->havePdbInfo = bPdbinfo;
    if (atoms->havePdbInfo)
    {
        snew(atoms->pdbinfo, natoms);
    }
    else
    {
        atoms->pdbinfo = nullptr;
    }
}

void gmx_pdbinfo_init_default(t_pdbinfo* pdbinfo)
{
    pdbinfo->type         = PdbRecordType::Atom;
    pdbinfo->atomnr       = 0;
    pdbinfo->altloc       = ' ';
    pdbinfo->atomnm[0]    = '\0';
    pdbinfo->occup        = 1.0;
    pdbinfo->bfac         = 0.0;
    pdbinfo->bAnisotropic = FALSE;
    std::fill(pdbinfo->uij, pdbinfo->uij + 6, 0.0);
}

t_atoms* copy_t_atoms(const t_atoms* src)
{
    t_atoms* dst = nullptr;

    snew(dst, 1);
    init_t_atoms(dst, src->nr, (nullptr != src->pdbinfo));
    if (nullptr != src->atomtype)
    {
        snew(dst->atomtype, src->nr);
    }
    if (nullptr != src->atomtypeB)
    {
        snew(dst->atomtypeB, src->nr);
    }
    for (int i = 0; (i < src->nr); i++)
    {
        dst->atom[i] = src->atom[i];
        if (nullptr != src->pdbinfo)
        {
            dst->pdbinfo[i] = src->pdbinfo[i];
        }
        if (nullptr != src->atomname)
        {
            dst->atomname[i] = src->atomname[i];
        }
        if (nullptr != src->atomtype)
        {
            dst->atomtype[i] = src->atomtype[i];
        }
        if (nullptr != src->atomtypeB)
        {
            dst->atomtypeB[i] = src->atomtypeB[i];
        }
    }
    dst->haveBState  = src->haveBState;
    dst->haveCharge  = src->haveCharge;
    dst->haveMass    = src->haveMass;
    dst->havePdbInfo = src->havePdbInfo;
    dst->haveType    = src->haveType;
    dst->nres        = src->nres;
    for (int i = 0; (i < src->nres); i++)
    {
        dst->resinfo[i] = src->resinfo[i];
    }
    return dst;
}

void t_atoms_set_resinfo(t_atoms*      atoms,
                         int           atom_ind,
                         t_symtab*     symtab,
                         const char*   resname,
                         int           resnr,
                         unsigned char ic,
                         int           chainnum,
                         char          chainid)
{
    t_resinfo* ri = &atoms->resinfo[atoms->atom[atom_ind].resind];
    ri->name      = put_symtab(symtab, resname);
    ri->rtp       = nullptr;
    ri->nr        = resnr;
    ri->ic        = ic;
    ri->chainnum  = chainnum;
    ri->chainid   = chainid;
}

static void pr_atom(FILE* fp, int indent, const char* title, const t_atom* atom, int n)
{
    if (available(fp, atom, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (int i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp,
                    "%s[%6d]={type=%3hu, typeB=%3hu, ptype=%8s, m=%12.5e, "
                    "q=%12.5e, mB=%12.5e, qB=%12.5e, resind=%5d, atomnumber=%3d}\n",
                    title,
                    i,
                    atom[i].type,
                    atom[i].typeB,
                    enumValueToString(atom[i].ptype),
                    atom[i].m,
                    atom[i].q,
                    atom[i].mB,
                    atom[i].qB,
                    atom[i].resind,
                    atom[i].atomnumber);
        }
    }
}

static void pr_strings2(FILE* fp, int indent, const char* title, char*** nm, char*** nmB, int n, gmx_bool bShowNumbers)
{
    if (available(fp, nm, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (int i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\",nameB=\"%s\"}\n", title, bShowNumbers ? i : -1, *(nm[i]), *(nmB[i]));
        }
    }
}

static void pr_resinfo(FILE* fp, int indent, const char* title, const t_resinfo* resinfo, int n, gmx_bool bShowNumbers)
{
    if (available(fp, resinfo, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (int i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp,
                    "%s[%d]={name=\"%s\", nr=%d, ic='%c'}\n",
                    title,
                    bShowNumbers ? i : -1,
                    *(resinfo[i].name),
                    resinfo[i].nr,
                    (resinfo[i].ic == '\0') ? ' ' : resinfo[i].ic);
        }
    }
}

void pr_atoms(FILE* fp, int indent, const char* title, const t_atoms* atoms, gmx_bool bShownumbers)
{
    if (available(fp, atoms, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_atom(fp, indent, "atom", atoms->atom, atoms->nr);
        pr_strings(fp, indent, "atom", atoms->atomname, atoms->nr, bShownumbers);
        pr_strings2(fp, indent, "type", atoms->atomtype, atoms->atomtypeB, atoms->nr, bShownumbers);
        pr_resinfo(fp, indent, "residue", atoms->resinfo, atoms->nres, bShownumbers);
    }
}

static void compareAtom(FILE* fp, int index, const t_atom* a1, const t_atom* a2, real relativeTolerance, real absoluteTolerance)
{
    if (a2)
    {
        cmp_us(fp, "atom.type", index, a1->type, a2->type);
        cmpEnum<ParticleType>(fp, "atom.ptype", a1->ptype, a2->ptype);
        cmp_int(fp, "atom.resind", index, a1->resind, a2->resind);
        cmp_int(fp, "atom.atomnumber", index, a1->atomnumber, a2->atomnumber);
        cmp_real(fp, "atom.m", index, a1->m, a2->m, relativeTolerance, absoluteTolerance);
        cmp_real(fp, "atom.q", index, a1->q, a2->q, relativeTolerance, absoluteTolerance);
        cmp_us(fp, "atom.typeB", index, a1->typeB, a2->typeB);
        cmp_real(fp, "atom.mB", index, a1->mB, a2->mB, relativeTolerance, absoluteTolerance);
        cmp_real(fp, "atom.qB", index, a1->qB, a2->qB, relativeTolerance, absoluteTolerance);
        cmp_str(fp, "elem", index, a1->elem, a2->elem);
    }
    else
    {
        cmp_us(fp, "atom.type", index, a1->type, a1->typeB);
        cmp_real(fp, "atom.m", index, a1->m, a1->mB, relativeTolerance, absoluteTolerance);
        cmp_real(fp, "atom.q", index, a1->q, a1->qB, relativeTolerance, absoluteTolerance);
    }
}

static void compareResinfo(FILE* fp, int residue, const t_resinfo& r1, const t_resinfo& r2)
{
    fprintf(fp, "comparing t_resinfo\n");
    cmp_str(fp, "name", residue, *r1.name, *r2.name);
    cmp_int(fp, "nr", residue, r1.nr, r2.nr);
    cmp_uc(fp, "ic", residue, r1.ic, r2.ic);
    cmp_int(fp, "chainnum", residue, r1.chainnum, r2.chainnum);
    cmp_uc(fp, "chainid", residue, r1.chainid, r2.chainid);
    if ((r1.rtp || r2.rtp) && (!r1.rtp || !r2.rtp))
    {
        fprintf(fp, "rtp info is present in topology %d but not in the other\n", r1.rtp ? 1 : 2);
    }
    if (r1.rtp && r2.rtp)
    {
        cmp_str(fp, "rtp", residue, *r1.rtp, *r2.rtp);
    }
}

static void comparePdbinfo(FILE*            fp,
                           int              pdb,
                           const t_pdbinfo& pdb1,
                           const t_pdbinfo& pdb2,
                           real             relativeTolerance,
                           real             absoluteTolerance)
{
    fprintf(fp, "comparing t_pdbinfo\n");
    cmpEnum<PdbRecordType>(fp, "type", pdb1.type, pdb2.type);
    cmp_int(fp, "atomnr", pdb, pdb1.atomnr, pdb2.atomnr);
    cmp_uc(fp, "altloc", pdb, pdb1.altloc, pdb2.altloc);
    cmp_str(fp, "atomnm", pdb, pdb1.atomnm, pdb2.atomnm);
    cmp_real(fp, "occup", pdb, pdb1.occup, pdb2.occup, relativeTolerance, absoluteTolerance);
    cmp_real(fp, "bfac", pdb, pdb1.bfac, pdb2.bfac, relativeTolerance, absoluteTolerance);
    cmp_bool(fp, "bAnistropic", pdb, pdb1.bAnisotropic, pdb2.bAnisotropic);
    for (int i = 0; i < 6; i++)
    {
        std::string buf = gmx::formatString("uij[%d]", i);
        cmp_int(fp, buf.c_str(), pdb, pdb1.uij[i], pdb2.uij[i]);
    }
}


void compareAtoms(FILE* fp, const t_atoms* a1, const t_atoms* a2, real relativeTolerance, real absoluteTolerance)
{
    fprintf(fp, "comparing atoms\n");

    if (a2)
    {
        cmp_int(fp, "atoms->nr", -1, a1->nr, a2->nr);
        cmp_int(fp, "atoms->nres", -1, a1->nres, a2->nres);
        cmp_bool(fp, "atoms->haveMass", -1, a1->haveMass, a2->haveMass);
        cmp_bool(fp, "atoms->haveCharge", -1, a1->haveCharge, a2->haveCharge);
        cmp_bool(fp, "atoms->haveType", -1, a1->haveType, a2->haveType);
        cmp_bool(fp, "atoms->haveBState", -1, a1->haveBState, a2->haveBState);
        cmp_bool(fp, "atoms->havePdbInfo", -1, a1->havePdbInfo, a2->havePdbInfo);
        for (int i = 0; i < std::min(a1->nr, a2->nr); i++)
        {
            compareAtom(fp, i, &(a1->atom[i]), &(a2->atom[i]), relativeTolerance, absoluteTolerance);
            if (a1->atomname && a2->atomname)
            {
                cmp_str(fp, "atomname", i, *a1->atomname[i], *a2->atomname[i]);
            }
            if (a1->havePdbInfo && a2->havePdbInfo)
            {
                comparePdbinfo(fp, i, a1->pdbinfo[i], a2->pdbinfo[i], relativeTolerance, absoluteTolerance);
            }
            if (a1->haveType && a2->haveType)
            {
                cmp_str(fp, "atomtype", i, *a1->atomtype[i], *a2->atomtype[i]);
            }
            if (a1->haveBState && a2->haveBState)
            {
                cmp_str(fp, "atomtypeB", i, *a1->atomtypeB[i], *a2->atomtypeB[i]);
            }
        }
        for (int i = 0; i < std::min(a1->nres, a2->nres); i++)
        {
            compareResinfo(fp, i, a1->resinfo[i], a2->resinfo[i]);
        }
    }
    else
    {
        for (int i = 0; (i < a1->nr); i++)
        {
            compareAtom(fp, i, &(a1->atom[i]), nullptr, relativeTolerance, absoluteTolerance);
        }
    }
}

void atomsSetMassesBasedOnNames(t_atoms* atoms, gmx_bool printMissingMasses)
{
    if (atoms->haveMass)
    {
        /* We could decide to anyhow assign then or generate a fatal error,
         * but it's probably most useful to keep the masses we have.
         */
        return;
    }

    int maxWarn = (printMissingMasses ? 10 : 0);
    int numWarn = 0;

    AtomProperties aps;

    bool haveMass = true;
    for (int i = 0; i < atoms->nr; i++)
    {
        if (!aps.setAtomProperty(epropMass,
                                 *atoms->resinfo[atoms->atom[i].resind].name,
                                 *atoms->atomname[i],
                                 &atoms->atom[i].m))
        {
            haveMass = false;

            if (numWarn < maxWarn)
            {
                fprintf(stderr,
                        "Can not find mass in database for atom %s in residue %d %s\n",
                        *atoms->atomname[i],
                        atoms->resinfo[atoms->atom[i].resind].nr,
                        *atoms->resinfo[atoms->atom[i].resind].name);
                numWarn++;
            }
            else
            {
                break;
            }
        }
    }
    atoms->haveMass = haveMass;
}
