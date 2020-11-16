/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/***************************************************************************
 *
 * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <algorithm>
#include <iostream>
#include <iomanip>

#include "starfile_label.h"
#include "starfile_error.h"

//This is needed for static memory allocation
std::map<EMDLabel, EMDLabelData> EMDL::data;
std::map<std::string, EMDLabel> EMDL::names;
std::map<std::string, std::string> EMDL::definitions;
StaticInitialization EMDL::initialization; //Just for initialization

void EMDL::addLabel(EMDLabel label, EMDLabelType type, std::string name, std::string definition)
{
    data[label] = EMDLabelData(type, name);
    names[name] = label;
    definitions[name] = definition;
}

void EMDL::addAltLabel(EMDLabel label, std::string name)
{
    names[name] = label;
}

void EMDL::printDefinitions(std::ostream& out)
{
    out << "+++ MetaDataLabel (EMDL) definitions: +++" <<std::endl;
    std::map<std::string, std::string>::const_iterator strIt;
    for (strIt = definitions.begin(); strIt != definitions.end(); strIt++)
    {
        out << std::setw(30) <<strIt->first;
        if (EMDL::isInt(names[strIt->first]))
            out << " (int)    ";
        else if (EMDL::isBool(names[strIt->first]))
            out << " (bool)   ";
        else if (EMDL::isDouble(names[strIt->first]))
            out << " (double) ";
        else if (EMDL::isString(names[strIt->first]))
            out << " (string) ";
        else
            REPORT_ERROR("EMDL::printDefinitions: unrecognised type");
        out << ": " << strIt->second <<std::endl;
    }
}



EMDLabel  EMDL::str2Label(const std::string &labelName)
{
    if (names.find(labelName) == names.end())
        return EMDL_UNDEFINED;
    return names[labelName];
}//close function str2Label

std::string  EMDL::label2Str(const EMDLabel &label)
{
    if (data.find(label) == data.end())
        return "";
    return data[label].str;
}//close function label2Str

bool EMDL::isInt(const EMDLabel &label)
{
    return (data[label].type == EMDL_INT);
}
bool EMDL::isBool(const EMDLabel &label)
{
    return (data[label].type == EMDL_BOOL);
}
bool EMDL::isString(const EMDLabel &label)
{
    return (data[label].type == EMDL_STRING);
}
bool EMDL::isDouble(const EMDLabel &label)
{
    return (data[label].type == EMDL_DOUBLE);
}
bool EMDL::isNumber(const EMDLabel &label)
{
    return (data[label].type == EMDL_DOUBLE || data[label].type == EMDL_INT);
}

bool EMDL::isValidLabel(const EMDLabel &label)
{
    return (label > EMDL_UNDEFINED && label < EMDL_LAST_LABEL);
}
bool EMDL::isValidLabel(const std::string &labelName)
{
    EMDLabel label = EMDL::str2Label(labelName);
    return EMDL::isValidLabel(label);
}

bool vectorContainsLabel(const std::vector<EMDLabel>& labelsVector, const EMDLabel label)
{
    std::vector<EMDLabel>::const_iterator location;
    location = std::find(labelsVector.begin(), labelsVector.end(), label);

    return (location != labelsVector.end());
}

/* in case we ever need to print out the enum values for verification
void
verify_enum()
{
  using namespace std;
  cout << "EMDL_UNDEFINED = " << EMDL_UNDEFINED << endl;
  cout << "EMDL_FIRST_LABEL = " << EMDL_FIRST_LABEL << endl;
  cout << "EMDL_OBJID = " << EMDL_OBJID << endl;
  cout << "EMDL_AREA_ID = " << EMDL_AREA_ID << endl;
  cout << "EMDL_AREA_NAME = " << EMDL_AREA_NAME << endl;
  cout << "EMDL_COMMENT = " << EMDL_COMMENT << endl;
  cout << "EMDL_BODY_MASK_NAME = " << EMDL_BODY_MASK_NAME << endl;
  cout << "EMDL_BODY_KEEP_FIXED = " << EMDL_BODY_KEEP_FIXED << endl;
  cout << "EMDL_BODY_REFERENCE_NAME = " << EMDL_BODY_REFERENCE_NAME << endl;
  cout << "EMDL_BODY_ROTATE_DIRECTION_X = " << EMDL_BODY_ROTATE_DIRECTION_X << endl;
  cout << "EMDL_BODY_ROTATE_DIRECTION_Y = " << EMDL_BODY_ROTATE_DIRECTION_Y << endl;
  cout << "EMDL_BODY_ROTATE_DIRECTION_Z = " << EMDL_BODY_ROTATE_DIRECTION_Z << endl;
  cout << "EMDL_BODY_ROTATE_RELATIVE_TO = " << EMDL_BODY_ROTATE_RELATIVE_TO << endl;
  cout << "EMDL_BODY_SIGMA_ANG = " << EMDL_BODY_SIGMA_ANG << endl;
  cout << "EMDL_BODY_SIGMA_OFFSET = " << EMDL_BODY_SIGMA_OFFSET << endl;
  cout << "EMDL_BODY_SIGMA_ROT = " << EMDL_BODY_SIGMA_ROT << endl;
  cout << "EMDL_BODY_SIGMA_TILT = " << EMDL_BODY_SIGMA_TILT << endl;
  cout << "EMDL_BODY_SIGMA_PSI = " << EMDL_BODY_SIGMA_PSI << endl;
  cout << "EMDL_BODY_STAR_FILE = " << EMDL_BODY_STAR_FILE << endl;
  cout << "EMDL_CTF_ASTIGMATISM = " << EMDL_CTF_ASTIGMATISM << endl;
  cout << "EMDL_CTF_BFACTOR = " << EMDL_CTF_BFACTOR << endl;
  cout << "EMDL_CTF_MAXRES = " << EMDL_CTF_MAXRES << endl;
  cout << "EMDL_CTF_VALIDATIONSCORE = " << EMDL_CTF_VALIDATIONSCORE << endl;
  cout << "EMDL_CTF_SCALEFACTOR = " << EMDL_CTF_SCALEFACTOR << endl;
  cout << "EMDL_CTF_SAMPLING_RATE = " << EMDL_CTF_SAMPLING_RATE << endl;
  cout << "EMDL_CTF_VOLTAGE = " << EMDL_CTF_VOLTAGE << endl;
  cout << "EMDL_CTF_DEFOCUSU = " << EMDL_CTF_DEFOCUSU << endl;
  cout << "EMDL_CTF_DEFOCUSV = " << EMDL_CTF_DEFOCUSV << endl;
  cout << "EMDL_CTF_DEFOCUS_ANGLE = " << EMDL_CTF_DEFOCUS_ANGLE << endl;
  cout << "EMDL_CTF_CS = " << EMDL_CTF_CS << endl;
  cout << "EMDL_CTF_CA = " << EMDL_CTF_CA << endl;
  cout << "EMDL_CTF_DETECTOR_PIXEL_SIZE = " << EMDL_CTF_DETECTOR_PIXEL_SIZE << endl;
  cout << "EMDL_CTF_ENERGY_LOSS = " << EMDL_CTF_ENERGY_LOSS << endl;
  cout << "EMDL_CTF_IMAGE = " << EMDL_CTF_IMAGE << endl;
  cout << "EMDL_CTF_LENS_STABILITY = " << EMDL_CTF_LENS_STABILITY << endl;
  cout << "EMDL_CTF_MAGNIFICATION = " << EMDL_CTF_MAGNIFICATION << endl;
  cout << "EMDL_CTF_PHASESHIFT = " << EMDL_CTF_PHASESHIFT << endl;
  cout << "EMDL_CTF_CONVERGENCE_CONE = " << EMDL_CTF_CONVERGENCE_CONE << endl;
  cout << "EMDL_CTF_LONGITUDINAL_DISPLACEMENT = " << EMDL_CTF_LONGITUDINAL_DISPLACEMENT << endl;
  cout << "EMDL_CTF_TRANSVERSAL_DISPLACEMENT = " << EMDL_CTF_TRANSVERSAL_DISPLACEMENT << endl;
  cout << "EMDL_CTF_Q0 = " << EMDL_CTF_Q0 << endl;
  cout << "EMDL_CTF_K = " << EMDL_CTF_K << endl;
  cout << "EMDL_CTF_VALUE = " << EMDL_CTF_VALUE << endl;
  cout << "EMDL_IMAGE_NAME = " << EMDL_IMAGE_NAME << endl;
  cout << "EMDL_IMAGE_ORI_NAME = " << EMDL_IMAGE_ORI_NAME << endl;
  cout << "EMDL_IMAGE_RECONSTRUCT_NAME = " << EMDL_IMAGE_RECONSTRUCT_NAME << endl;
  cout << "EMDL_IMAGE_ID = " << EMDL_IMAGE_ID << endl;
  cout << "EMDL_IMAGE_DATATYPE = " << EMDL_IMAGE_DATATYPE << endl;
  cout << "EMDL_IMAGE_DIMENSIONALITY = " << EMDL_IMAGE_DIMENSIONALITY << endl;
  cout << "EMDL_IMAGE_BEAMTILT_X = " << EMDL_IMAGE_BEAMTILT_X << endl;
  cout << "EMDL_IMAGE_BEAMTILT_Y = " << EMDL_IMAGE_BEAMTILT_Y << endl;
  cout << "EMDL_IMAGE_COORD_X = " << EMDL_IMAGE_COORD_X << endl;
  cout << "EMDL_IMAGE_COORD_Y = " << EMDL_IMAGE_COORD_Y << endl;
  cout << "EMDL_IMAGE_COORD_Z = " << EMDL_IMAGE_COORD_Z << endl;
  cout << "EMDL_IMAGE_FRAME_NR = " << EMDL_IMAGE_FRAME_NR << endl;
  cout << "EMDL_IMAGE_MAGNIFICATION_CORRECTION = " << EMDL_IMAGE_MAGNIFICATION_CORRECTION << endl;
  cout << "EMDL_IMAGE_NORM_CORRECTION = " << EMDL_IMAGE_NORM_CORRECTION << endl;
  cout << "EMDL_IMAGE_SAMPLINGRATE = " << EMDL_IMAGE_SAMPLINGRATE << endl;
  cout << "EMDL_IMAGE_SAMPLINGRATE_X = " << EMDL_IMAGE_SAMPLINGRATE_X << endl;
  cout << "EMDL_IMAGE_SAMPLINGRATE_Y = " << EMDL_IMAGE_SAMPLINGRATE_Y << endl;
  cout << "EMDL_IMAGE_SAMPLINGRATE_Z = " << EMDL_IMAGE_SAMPLINGRATE_Z << endl;
  cout << "EMDL_IMAGE_SIZE = " << EMDL_IMAGE_SIZE << endl;
  cout << "EMDL_IMAGE_SIZE_X = " << EMDL_IMAGE_SIZE_X << endl;
  cout << "EMDL_IMAGE_SIZE_Y = " << EMDL_IMAGE_SIZE_Y << endl;
  cout << "EMDL_IMAGE_SIZE_Z = " << EMDL_IMAGE_SIZE_Z << endl;
  cout << "EMDL_IMAGE_STATS_MIN = " << EMDL_IMAGE_STATS_MIN << endl;
  cout << "EMDL_IMAGE_STATS_MAX = " << EMDL_IMAGE_STATS_MAX << endl;
  cout << "EMDL_IMAGE_STATS_AVG = " << EMDL_IMAGE_STATS_AVG << endl;
  cout << "EMDL_IMAGE_STATS_STDDEV = " << EMDL_IMAGE_STATS_STDDEV << endl;
  cout << "EMDL_IMAGE_STATS_SKEW = " << EMDL_IMAGE_STATS_SKEW << endl;
  cout << "EMDL_IMAGE_STATS_KURT = " << EMDL_IMAGE_STATS_KURT << endl;
  cout << "EMDL_IMAGE_WEIGHT = " << EMDL_IMAGE_WEIGHT << endl;
  cout << "EMDL_MATRIX_1_1 = " << EMDL_MATRIX_1_1 << endl;
  cout << "EMDL_MATRIX_1_2 = " << EMDL_MATRIX_1_2 << endl;
  cout << "EMDL_MATRIX_1_3 = " << EMDL_MATRIX_1_3 << endl;
  cout << "EMDL_MATRIX_2_1 = " << EMDL_MATRIX_2_1 << endl;
  cout << "EMDL_MATRIX_2_2 = " << EMDL_MATRIX_2_2 << endl;
  cout << "EMDL_MATRIX_2_3 = " << EMDL_MATRIX_2_3 << endl;
  cout << "EMDL_MATRIX_3_1 = " << EMDL_MATRIX_3_1 << endl;
  cout << "EMDL_MATRIX_3_2 = " << EMDL_MATRIX_3_2 << endl;
  cout << "EMDL_MATRIX_3_3 = " << EMDL_MATRIX_3_3 << endl;
  cout << "EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL = " << EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL << endl;
  cout << "EMDL_MICROGRAPH_ACCUM_MOTION_EARLY = " << EMDL_MICROGRAPH_ACCUM_MOTION_EARLY << endl;
  cout << "EMDL_MICROGRAPH_ACCUM_MOTION_LATE = " << EMDL_MICROGRAPH_ACCUM_MOTION_LATE << endl;
  cout << "EMDL_MICROGRAPH_ID = " << EMDL_MICROGRAPH_ID << endl;
  cout << "EMDL_MICROGRAPH_NAME = " << EMDL_MICROGRAPH_NAME << endl;
  cout << "EMDL_MICROGRAPH_GAIN_NAME = " << EMDL_MICROGRAPH_GAIN_NAME << endl;
  cout << "EMDL_MICROGRAPH_DEFECT_FILE = " << EMDL_MICROGRAPH_DEFECT_FILE << endl;
  cout << "EMDL_MICROGRAPH_NAME_WODOSE = " << EMDL_MICROGRAPH_NAME_WODOSE << endl;
  cout << "EMDL_MICROGRAPH_MOVIE_NAME = " << EMDL_MICROGRAPH_MOVIE_NAME << endl;
  cout << "EMDL_MICROGRAPH_METADATA_NAME = " << EMDL_MICROGRAPH_METADATA_NAME << endl;
  cout << "EMDL_MICROGRAPH_TILT_ANGLE = " << EMDL_MICROGRAPH_TILT_ANGLE << endl;
  cout << "EMDL_MICROGRAPH_TILT_AXIS_DIRECTION = " << EMDL_MICROGRAPH_TILT_AXIS_DIRECTION << endl;
  cout << "EMDL_MICROGRAPH_TILT_AXIS_OUTOFPLANE = " << EMDL_MICROGRAPH_TILT_AXIS_OUTOFPLANE << endl;
  cout << "EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE = " << EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE << endl;
  cout << "EMDL_MICROGRAPH_PRE_EXPOSURE = " << EMDL_MICROGRAPH_PRE_EXPOSURE << endl;
  cout << "EMDL_MICROGRAPH_DOSE_RATE = " << EMDL_MICROGRAPH_DOSE_RATE << endl;
  cout << "EMDL_MICROGRAPH_BINNING = " << EMDL_MICROGRAPH_BINNING << endl;
  cout << "EMDL_MICROGRAPH_FRAME_NUMBER = " << EMDL_MICROGRAPH_FRAME_NUMBER << endl;
  cout << "EMDL_MICROGRAPH_MOTION_MODEL_VERSION = " << EMDL_MICROGRAPH_MOTION_MODEL_VERSION << endl;
  cout << "EMDL_MICROGRAPH_START_FRAME = " << EMDL_MICROGRAPH_START_FRAME << endl;
  cout << "EMDL_MICROGRAPH_END_FRAME = " << EMDL_MICROGRAPH_END_FRAME << endl;
  cout << "EMDL_MICROGRAPH_SHIFT_X = " << EMDL_MICROGRAPH_SHIFT_X << endl;
  cout << "EMDL_MICROGRAPH_SHIFT_Y = " << EMDL_MICROGRAPH_SHIFT_Y << endl;
  cout << "EMDL_MICROGRAPH_MOTION_COEFFS_IDX = " << EMDL_MICROGRAPH_MOTION_COEFFS_IDX << endl;
  cout << "EMDL_MICROGRAPH_MOTION_COEFF = " << EMDL_MICROGRAPH_MOTION_COEFF << endl;
  cout << "EMDL_MASK_NAME = " << EMDL_MASK_NAME << endl;
  cout << "EMDL_MLMODEL_ACCURACY_ROT = " << EMDL_MLMODEL_ACCURACY_ROT << endl;
  cout << "EMDL_MLMODEL_ACCURACY_TRANS = " << EMDL_MLMODEL_ACCURACY_TRANS << endl;
  cout << "EMDL_MLMODEL_AVE_PMAX = " << EMDL_MLMODEL_AVE_PMAX << endl;
  cout << "EMDL_MLMODEL_CURRENT_RESOLUTION = " << EMDL_MLMODEL_CURRENT_RESOLUTION << endl;
  cout << "EMDL_MLMODEL_CURRENT_SIZE = " << EMDL_MLMODEL_CURRENT_SIZE << endl;
  cout << "EMDL_MLMODEL_DATA_VS_PRIOR_REF = " << EMDL_MLMODEL_DATA_VS_PRIOR_REF << endl;
  cout << "EMDL_MLMODEL_DIMENSIONALITY = " << EMDL_MLMODEL_DIMENSIONALITY << endl;
  cout << "EMDL_MLMODEL_DIMENSIONALITY_DATA = " << EMDL_MLMODEL_DIMENSIONALITY_DATA << endl;
  cout << "EMDL_MLMODEL_DIFF2_HALVES_REF = " << EMDL_MLMODEL_DIFF2_HALVES_REF << endl;
  cout << "EMDL_MLMODEL_ESTIM_RESOL_REF = " << EMDL_MLMODEL_ESTIM_RESOL_REF << endl;
  cout << "EMDL_MLMODEL_FOURIER_COVERAGE_REF = " << EMDL_MLMODEL_FOURIER_COVERAGE_REF << endl;
  cout << "EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF = " << EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF << endl;
  cout << "EMDL_MLMODEL_FSC_HALVES_REF = " << EMDL_MLMODEL_FSC_HALVES_REF << endl;
  cout << "EMDL_MLMODEL_GROUP_NAME = " << EMDL_MLMODEL_GROUP_NAME << endl;
  cout << "EMDL_MLMODEL_GROUP_NO = " << EMDL_MLMODEL_GROUP_NO << endl;
  cout << "EMDL_MLMODEL_GROUP_NR_PARTICLES = " << EMDL_MLMODEL_GROUP_NR_PARTICLES << endl;
  cout << "EMDL_MLMODEL_GROUP_SCALE_CORRECTION = " << EMDL_MLMODEL_GROUP_SCALE_CORRECTION << endl;
  cout << "EMDL_MLMODEL_HELICAL_NR_ASU = " << EMDL_MLMODEL_HELICAL_NR_ASU << endl;
  cout << "EMDL_MLMODEL_HELICAL_TWIST = " << EMDL_MLMODEL_HELICAL_TWIST << endl;
  cout << "EMDL_MLMODEL_HELICAL_TWIST_MIN = " << EMDL_MLMODEL_HELICAL_TWIST_MIN << endl;
  cout << "EMDL_MLMODEL_HELICAL_TWIST_MAX = " << EMDL_MLMODEL_HELICAL_TWIST_MAX << endl;
  cout << "EMDL_MLMODEL_HELICAL_TWIST_INITIAL_STEP = " << EMDL_MLMODEL_HELICAL_TWIST_INITIAL_STEP << endl;
  cout << "EMDL_MLMODEL_HELICAL_RISE = " << EMDL_MLMODEL_HELICAL_RISE << endl;
  cout << "EMDL_MLMODEL_HELICAL_RISE_MIN = " << EMDL_MLMODEL_HELICAL_RISE_MIN << endl;
  cout << "EMDL_MLMODEL_HELICAL_RISE_MAX = " << EMDL_MLMODEL_HELICAL_RISE_MAX << endl;
  cout << "EMDL_MLMODEL_HELICAL_RISE_INITIAL_STEP = " << EMDL_MLMODEL_HELICAL_RISE_INITIAL_STEP << endl;
  cout << "EMDL_MLMODEL_IS_HELIX = " << EMDL_MLMODEL_IS_HELIX << endl;
  cout << "EMDL_MLMODEL_INTERPOLATOR = " << EMDL_MLMODEL_INTERPOLATOR << endl;
  cout << "EMDL_MLMODEL_LL = " << EMDL_MLMODEL_LL << endl;
  cout << "EMDL_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION = " << EMDL_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION << endl;
  cout << "EMDL_MLMODEL_NORM_CORRECTION_AVG = " << EMDL_MLMODEL_NORM_CORRECTION_AVG << endl;
  cout << "EMDL_MLMODEL_NR_BODIES = " << EMDL_MLMODEL_NR_BODIES << endl;
  cout << "EMDL_MLMODEL_NR_CLASSES = " << EMDL_MLMODEL_NR_CLASSES << endl;
  cout << "EMDL_MLMODEL_NR_GROUPS = " << EMDL_MLMODEL_NR_GROUPS << endl;
  cout << "EMDL_MLMODEL_ORIGINAL_SIZE = " << EMDL_MLMODEL_ORIGINAL_SIZE << endl;
  cout << "EMDL_MLMODEL_ORIENTABILITY_CONTRIBUTION = " << EMDL_MLMODEL_ORIENTABILITY_CONTRIBUTION << endl;
  cout << "EMDL_MLMODEL_PADDING_FACTOR = " << EMDL_MLMODEL_PADDING_FACTOR << endl;
  cout << "EMDL_MLMODEL_PDF_CLASS = " << EMDL_MLMODEL_PDF_CLASS << endl;
  cout << "EMDL_MLMODEL_PRIOR_OFFX_CLASS = " << EMDL_MLMODEL_PRIOR_OFFX_CLASS << endl;
  cout << "EMDL_MLMODEL_PRIOR_OFFY_CLASS = " << EMDL_MLMODEL_PRIOR_OFFY_CLASS << endl;
  cout << "EMDL_MLMODEL_PDF_ORIENT = " << EMDL_MLMODEL_PDF_ORIENT << endl;
  cout << "EMDL_MLMODEL_PIXEL_SIZE = " << EMDL_MLMODEL_PIXEL_SIZE << endl;
  cout << "EMDL_MLMODEL_POWER_REF = " << EMDL_MLMODEL_POWER_REF << endl;
  cout << "EMDL_MLMODEL_PRIOR_MODE = " << EMDL_MLMODEL_PRIOR_MODE << endl;
  cout << "EMDL_MLMODEL_SIGMA_OFFSET = " << EMDL_MLMODEL_SIGMA_OFFSET << endl;
  cout << "EMDL_MLMODEL_SIGMA_ROT = " << EMDL_MLMODEL_SIGMA_ROT << endl;
  cout << "EMDL_MLMODEL_SIGMA_TILT = " << EMDL_MLMODEL_SIGMA_TILT << endl;
  cout << "EMDL_MLMODEL_SIGMA_PSI = " << EMDL_MLMODEL_SIGMA_PSI << endl;
  cout << "EMDL_MLMODEL_REF_IMAGE = " << EMDL_MLMODEL_REF_IMAGE << endl;
  cout << "EMDL_MLMODEL_SGD_GRADIENT_IMAGE = " << EMDL_MLMODEL_SGD_GRADIENT_IMAGE << endl;
  cout << "EMDL_MLMODEL_SIGMA2_NOISE = " << EMDL_MLMODEL_SIGMA2_NOISE << endl;
  cout << "EMDL_MLMODEL_SIGMA2_REF = " << EMDL_MLMODEL_SIGMA2_REF << endl;
  cout << "EMDL_MLMODEL_SSNR_REF = " << EMDL_MLMODEL_SSNR_REF << endl;
  cout << "EMDL_MLMODEL_TAU2_FUDGE_FACTOR = " << EMDL_MLMODEL_TAU2_FUDGE_FACTOR << endl;
  cout << "EMDL_MLMODEL_TAU2_REF = " << EMDL_MLMODEL_TAU2_REF << endl;
  cout << "EMDL_OPTIMISER_ACCURACY_ROT = " << EMDL_OPTIMISER_ACCURACY_ROT << endl;
  cout << "EMDL_OPTIMISER_ACCURACY_TRANS = " << EMDL_OPTIMISER_ACCURACY_TRANS << endl;
  cout << "EMDL_OPTIMISER_ADAPTIVE_FRACTION = " << EMDL_OPTIMISER_ADAPTIVE_FRACTION << endl;
  cout << "EMDL_OPTIMISER_ADAPTIVE_OVERSAMPLING = " << EMDL_OPTIMISER_ADAPTIVE_OVERSAMPLING << endl;
  cout << "EMDL_OPTIMISER_AUTO_LOCAL_HP_ORDER = " << EMDL_OPTIMISER_AUTO_LOCAL_HP_ORDER << endl;
  cout << "EMDL_OPTIMISER_AVAILABLE_MEMORY = " << EMDL_OPTIMISER_AVAILABLE_MEMORY << endl;
  cout << "EMDL_OPTIMISER_BEST_RESOL_THUS_FAR = " << EMDL_OPTIMISER_BEST_RESOL_THUS_FAR << endl;
  cout << "EMDL_OPTIMISER_CHANGES_OPTIMAL_OFFSETS = " << EMDL_OPTIMISER_CHANGES_OPTIMAL_OFFSETS << endl;
  cout << "EMDL_OPTIMISER_CHANGES_OPTIMAL_ORIENTS = " << EMDL_OPTIMISER_CHANGES_OPTIMAL_ORIENTS << endl;
  cout << "EMDL_OPTIMISER_CHANGES_OPTIMAL_CLASSES = " << EMDL_OPTIMISER_CHANGES_OPTIMAL_CLASSES << endl;
  cout << "EMDL_OPTIMISER_COARSE_SIZE = " << EMDL_OPTIMISER_COARSE_SIZE << endl;
  cout << "EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED = " << EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED << endl;
  cout << "EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED = " << EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED << endl;
  cout << "EMDL_OPTIMISER_DATA_STARFILE = " << EMDL_OPTIMISER_DATA_STARFILE << endl;
  cout << "EMDL_OPTIMISER_DO_AUTO_REFINE = " << EMDL_OPTIMISER_DO_AUTO_REFINE << endl;
  cout << "EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES = " << EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES << endl;
  cout << "EMDL_OPTIMISER_DO_CORRECT_CTF = " << EMDL_OPTIMISER_DO_CORRECT_CTF << endl;
  cout << "EMDL_OPTIMISER_DO_CORRECT_MAGNIFICATION = " << EMDL_OPTIMISER_DO_CORRECT_MAGNIFICATION << endl;
  cout << "EMDL_OPTIMISER_DO_CORRECT_NORM = " << EMDL_OPTIMISER_DO_CORRECT_NORM << endl;
  cout << "EMDL_OPTIMISER_DO_CORRECT_SCALE = " << EMDL_OPTIMISER_DO_CORRECT_SCALE << endl;
  cout << "EMDL_OPTIMISER_DO_REALIGN_MOVIES = " << EMDL_OPTIMISER_DO_REALIGN_MOVIES << endl;
  cout << "EMDL_OPTIMISER_DO_MAP = " << EMDL_OPTIMISER_DO_MAP << endl;
  cout << "EMDL_OPTIMISER_DO_SGD = " << EMDL_OPTIMISER_DO_SGD << endl;
  cout << "EMDL_OPTIMISER_FAST_SUBSETS = " << EMDL_OPTIMISER_FAST_SUBSETS << endl;
  cout << "EMDL_OPTIMISER_SGD_INI_ITER = " << EMDL_OPTIMISER_SGD_INI_ITER << endl;
  cout << "EMDL_OPTIMISER_SGD_FIN_ITER = " << EMDL_OPTIMISER_SGD_FIN_ITER << endl;
  cout << "EMDL_OPTIMISER_SGD_INBETWEEN_ITER = " << EMDL_OPTIMISER_SGD_INBETWEEN_ITER << endl;
  cout << "EMDL_OPTIMISER_SGD_INI_RESOL = " << EMDL_OPTIMISER_SGD_INI_RESOL << endl;
  cout << "EMDL_OPTIMISER_SGD_FIN_RESOL = " << EMDL_OPTIMISER_SGD_FIN_RESOL << endl;
  cout << "EMDL_OPTIMISER_SGD_INI_SUBSET_SIZE = " << EMDL_OPTIMISER_SGD_INI_SUBSET_SIZE << endl;
  cout << "EMDL_OPTIMISER_SGD_FIN_SUBSET_SIZE = " << EMDL_OPTIMISER_SGD_FIN_SUBSET_SIZE << endl;
  cout << "EMDL_OPTIMISER_SGD_MU = " << EMDL_OPTIMISER_SGD_MU << endl;
  cout << "EMDL_OPTIMISER_SGD_SIGMA2FUDGE_INI = " << EMDL_OPTIMISER_SGD_SIGMA2FUDGE_INI << endl;
  cout << "EMDL_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE = " << EMDL_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE << endl;
  cout << "EMDL_OPTIMISER_SGD_SKIP_ANNNEAL = " << EMDL_OPTIMISER_SGD_SKIP_ANNNEAL << endl;
  cout << "EMDL_OPTIMISER_SGD_SUBSET_SIZE = " << EMDL_OPTIMISER_SGD_SUBSET_SIZE << endl;
  cout << "EMDL_OPTIMISER_SGD_WRITE_EVERY_SUBSET = " << EMDL_OPTIMISER_SGD_WRITE_EVERY_SUBSET << endl;
  cout << "EMDL_OPTIMISER_SGD_MAX_SUBSETS = " << EMDL_OPTIMISER_SGD_MAX_SUBSETS << endl;
  cout << "EMDL_OPTIMISER_SGD_STEPSIZE = " << EMDL_OPTIMISER_SGD_STEPSIZE << endl;
  cout << "EMDL_OPTIMISER_DO_SOLVENT_FLATTEN = " << EMDL_OPTIMISER_DO_SOLVENT_FLATTEN << endl;
  cout << "EMDL_OPTIMISER_DO_SOLVENT_FSC = " << EMDL_OPTIMISER_DO_SOLVENT_FSC << endl;
  cout << "EMDL_OPTIMISER_DO_SKIP_ALIGN = " << EMDL_OPTIMISER_DO_SKIP_ALIGN << endl;
  cout << "EMDL_OPTIMISER_DO_SKIP_ROTATE = " << EMDL_OPTIMISER_DO_SKIP_ROTATE << endl;
  cout << "EMDL_OPTIMISER_DO_SPLIT_RANDOM_HALVES = " << EMDL_OPTIMISER_DO_SPLIT_RANDOM_HALVES << endl;
  cout << "EMDL_OPTIMISER_DO_ZERO_MASK = " << EMDL_OPTIMISER_DO_ZERO_MASK << endl;
  cout << "EMDL_OPTIMISER_FIX_SIGMA_NOISE = " << EMDL_OPTIMISER_FIX_SIGMA_NOISE << endl;
  cout << "EMDL_OPTIMISER_FIX_SIGMA_OFFSET = " << EMDL_OPTIMISER_FIX_SIGMA_OFFSET << endl;
  cout << "EMDL_OPTIMISER_FIX_TAU = " << EMDL_OPTIMISER_FIX_TAU << endl;
  cout << "EMDL_OPTIMISER_HAS_CONVERGED = " << EMDL_OPTIMISER_HAS_CONVERGED << endl;
  cout << "EMDL_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT = " << EMDL_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT << endl;
  cout << "EMDL_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO = " << EMDL_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO << endl;
  cout << "EMDL_OPTIMISER_DO_HELICAL_REFINE = " << EMDL_OPTIMISER_DO_HELICAL_REFINE << endl;
  cout << "EMDL_OPTIMISER_IGNORE_HELICAL_SYMMETRY = " << EMDL_OPTIMISER_IGNORE_HELICAL_SYMMETRY << endl;
  cout << "EMDL_OPTIMISER_HELICAL_TWIST_INITIAL = " << EMDL_OPTIMISER_HELICAL_TWIST_INITIAL << endl;
  cout << "EMDL_OPTIMISER_HELICAL_RISE_INITIAL = " << EMDL_OPTIMISER_HELICAL_RISE_INITIAL << endl;
  cout << "EMDL_OPTIMISER_HELICAL_Z_PERCENTAGE = " << EMDL_OPTIMISER_HELICAL_Z_PERCENTAGE << endl;
  cout << "EMDL_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER = " << EMDL_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER << endl;
  cout << "EMDL_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER = " << EMDL_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER << endl;
  cout << "EMDL_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT = " << EMDL_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT << endl;
  cout << "EMDL_OPTIMISER_HELICAL_SIGMA_DISTANCE = " << EMDL_OPTIMISER_HELICAL_SIGMA_DISTANCE << endl;
  cout << "EMDL_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED = " << EMDL_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED << endl;
  cout << "EMDL_OPTIMISER_HIGHRES_LIMIT_EXP = " << EMDL_OPTIMISER_HIGHRES_LIMIT_EXP << endl;
  cout << "EMDL_OPTIMISER_HIGHRES_LIMIT_SGD = " << EMDL_OPTIMISER_HIGHRES_LIMIT_SGD << endl;
  cout << "EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK = " << EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK << endl;
  cout << "EMDL_OPTIMISER_INCR_SIZE = " << EMDL_OPTIMISER_INCR_SIZE << endl;
  cout << "EMDL_OPTIMISER_ITERATION_NO = " << EMDL_OPTIMISER_ITERATION_NO << endl;
  cout << "EMDL_OPTIMISER_LOCAL_SYMMETRY_FILENAME = " << EMDL_OPTIMISER_LOCAL_SYMMETRY_FILENAME << endl;
  cout << "EMDL_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES = " << EMDL_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES << endl;
  cout << "EMDL_OPTIMISER_MAGNIFICATION_RANGE = " << EMDL_OPTIMISER_MAGNIFICATION_RANGE << endl;
  cout << "EMDL_OPTIMISER_MAGNIFICATION_STEP = " << EMDL_OPTIMISER_MAGNIFICATION_STEP << endl;
  cout << "EMDL_OPTIMISER_MAX_COARSE_SIZE = " << EMDL_OPTIMISER_MAX_COARSE_SIZE << endl;
  cout << "EMDL_OPTIMISER_MAX_NR_POOL = " << EMDL_OPTIMISER_MAX_NR_POOL << endl;
  cout << "EMDL_OPTIMISER_MODEL_STARFILE = " << EMDL_OPTIMISER_MODEL_STARFILE << endl;
  cout << "EMDL_OPTIMISER_MODEL_STARFILE2 = " << EMDL_OPTIMISER_MODEL_STARFILE2 << endl;
  cout << "EMDL_OPTIMISER_NR_ITERATIONS = " << EMDL_OPTIMISER_NR_ITERATIONS << endl;
  cout << "EMDL_OPTIMISER_NR_ITER_WO_RESOL_GAIN = " << EMDL_OPTIMISER_NR_ITER_WO_RESOL_GAIN << endl;
  cout << "EMDL_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES = " << EMDL_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES << endl;
  cout << "EMDL_OPTIMISER_OUTPUT_ROOTNAME = " << EMDL_OPTIMISER_OUTPUT_ROOTNAME << endl;
  cout << "EMDL_OPTIMISER_PARTICLE_DIAMETER = " << EMDL_OPTIMISER_PARTICLE_DIAMETER << endl;
  cout << "EMDL_OPTIMISER_RADIUS_MASK_3D_MAP = " << EMDL_OPTIMISER_RADIUS_MASK_3D_MAP << endl;
  cout << "EMDL_OPTIMISER_RADIUS_MASK_EXP_PARTICLES = " << EMDL_OPTIMISER_RADIUS_MASK_EXP_PARTICLES << endl;
  cout << "EMDL_OPTIMISER_RANDOM_SEED = " << EMDL_OPTIMISER_RANDOM_SEED << endl;
  cout << "EMDL_OPTIMISER_REFS_ARE_CTF_CORRECTED = " << EMDL_OPTIMISER_REFS_ARE_CTF_CORRECTED << endl;
  cout << "EMDL_OPTIMISER_SAMPLING_STARFILE = " << EMDL_OPTIMISER_SAMPLING_STARFILE << endl;
  cout << "EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES = " << EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES << endl;
  cout << "EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS = " << EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS << endl;
  cout << "EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS = " << EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS << endl;
  cout << "EMDL_OPTIMISER_SOLVENT_MASK_NAME = " << EMDL_OPTIMISER_SOLVENT_MASK_NAME << endl;
  cout << "EMDL_OPTIMISER_SOLVENT_MASK2_NAME = " << EMDL_OPTIMISER_SOLVENT_MASK2_NAME << endl;
  cout << "EMDL_OPTIMISER_TAU_SPECTRUM_NAME = " << EMDL_OPTIMISER_TAU_SPECTRUM_NAME << endl;
  cout << "EMDL_OPTIMISER_USE_TOO_COARSE_SAMPLING = " << EMDL_OPTIMISER_USE_TOO_COARSE_SAMPLING << endl;
  cout << "EMDL_OPTIMISER_WIDTH_MASK_EDGE = " << EMDL_OPTIMISER_WIDTH_MASK_EDGE << endl;
  cout << "EMDL_ORIENT_FLIP = " << EMDL_ORIENT_FLIP << endl;
  cout << "EMDL_ORIENT_ID = " << EMDL_ORIENT_ID << endl;
  cout << "EMDL_ORIENT_ORIGIN_X = " << EMDL_ORIENT_ORIGIN_X << endl;
  cout << "EMDL_ORIENT_ORIGIN_X_PRIOR = " << EMDL_ORIENT_ORIGIN_X_PRIOR << endl;
  cout << "EMDL_ORIENT_ORIGIN_Y = " << EMDL_ORIENT_ORIGIN_Y << endl;
  cout << "EMDL_ORIENT_ORIGIN_Y_PRIOR = " << EMDL_ORIENT_ORIGIN_Y_PRIOR << endl;
  cout << "EMDL_ORIENT_ORIGIN_Z = " << EMDL_ORIENT_ORIGIN_Z << endl;
  cout << "EMDL_ORIENT_ORIGIN_Z_PRIOR = " << EMDL_ORIENT_ORIGIN_Z_PRIOR << endl;
  cout << "EMDL_ORIENT_ROT = " << EMDL_ORIENT_ROT << endl;
  cout << "EMDL_ORIENT_ROT_PRIOR = " << EMDL_ORIENT_ROT_PRIOR << endl;
  cout << "EMDL_ORIENT_TILT = " << EMDL_ORIENT_TILT << endl;
  cout << "EMDL_ORIENT_TILT_PRIOR = " << EMDL_ORIENT_TILT_PRIOR << endl;
  cout << "EMDL_ORIENT_PSI = " << EMDL_ORIENT_PSI << endl;
  cout << "EMDL_ORIENT_PSI_PRIOR = " << EMDL_ORIENT_PSI_PRIOR << endl;
  cout << "EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO = " << EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO << endl;
  cout << "EMDL_PARTICLE_AUTOPICK_FOM = " << EMDL_PARTICLE_AUTOPICK_FOM << endl;
  cout << "EMDL_PARTICLE_HELICAL_TUBE_ID = " << EMDL_PARTICLE_HELICAL_TUBE_ID << endl;
  cout << "EMDL_PARTICLE_HELICAL_TUBE_PITCH = " << EMDL_PARTICLE_HELICAL_TUBE_PITCH << endl;
  cout << "EMDL_PARTICLE_HELICAL_TRACK_LENGTH = " << EMDL_PARTICLE_HELICAL_TRACK_LENGTH << endl;
  cout << "EMDL_PARTICLE_CLASS = " << EMDL_PARTICLE_CLASS << endl;
  cout << "EMDL_PARTICLE_DLL = " << EMDL_PARTICLE_DLL << endl;
  cout << "EMDL_PARTICLE_ID = " << EMDL_PARTICLE_ID << endl;
  cout << "EMDL_PARTICLE_FOM = " << EMDL_PARTICLE_FOM << endl;
  cout << "EMDL_PARTICLE_KL_DIVERGENCE = " << EMDL_PARTICLE_KL_DIVERGENCE << endl;
  cout << "EMDL_PARTICLE_RANDOM_SUBSET = " << EMDL_PARTICLE_RANDOM_SUBSET << endl;
  cout << "EMDL_PARTICLE_BEAM_TILT_CLASS = " << EMDL_PARTICLE_BEAM_TILT_CLASS << endl;
  cout << "EMDL_PARTICLE_NAME = " << EMDL_PARTICLE_NAME << endl;
  cout << "EMDL_PARTICLE_ORI_NAME = " << EMDL_PARTICLE_ORI_NAME << endl;
  cout << "EMDL_PARTICLE_NR_SIGNIFICANT_SAMPLES = " << EMDL_PARTICLE_NR_SIGNIFICANT_SAMPLES << endl;
  cout << "EMDL_PARTICLE_NR_FRAMES = " << EMDL_PARTICLE_NR_FRAMES << endl;
  cout << "EMDL_PARTICLE_NR_FRAMES_AVG = " << EMDL_PARTICLE_NR_FRAMES_AVG << endl;
  cout << "EMDL_PARTICLE_MOVIE_RUNNING_AVG = " << EMDL_PARTICLE_MOVIE_RUNNING_AVG << endl;
  cout << "EMDL_PARTICLE_PMAX = " << EMDL_PARTICLE_PMAX << endl;
  cout << "EMDL_PARTICLE_NUMBER = " << EMDL_PARTICLE_NUMBER << endl;
  cout << "EMDL_PIPELINE_JOB_COUNTER = " << EMDL_PIPELINE_JOB_COUNTER << endl;
  cout << "EMDL_PIPELINE_NODE_NAME = " << EMDL_PIPELINE_NODE_NAME << endl;
  cout << "EMDL_PIPELINE_NODE_TYPE = " << EMDL_PIPELINE_NODE_TYPE << endl;
  cout << "EMDL_PIPELINE_PROCESS_ALIAS = " << EMDL_PIPELINE_PROCESS_ALIAS << endl;
  cout << "EMDL_PIPELINE_PROCESS_NAME = " << EMDL_PIPELINE_PROCESS_NAME << endl;
  cout << "EMDL_PIPELINE_PROCESS_TYPE = " << EMDL_PIPELINE_PROCESS_TYPE << endl;
  cout << "EMDL_PIPELINE_PROCESS_STATUS = " << EMDL_PIPELINE_PROCESS_STATUS << endl;
  cout << "EMDL_PIPELINE_EDGE_FROM = " << EMDL_PIPELINE_EDGE_FROM << endl;
  cout << "EMDL_PIPELINE_EDGE_TO = " << EMDL_PIPELINE_EDGE_TO << endl;
  cout << "EMDL_PIPELINE_EDGE_PROCESS = " << EMDL_PIPELINE_EDGE_PROCESS << endl;
  cout << "EMDL_POSTPROCESS_BFACTOR = " << EMDL_POSTPROCESS_BFACTOR << endl;
  cout << "EMDL_POSTPROCESS_FINAL_RESOLUTION = " << EMDL_POSTPROCESS_FINAL_RESOLUTION << endl;
  cout << "EMDL_POSTPROCESS_FSC_GENERAL = " << EMDL_POSTPROCESS_FSC_GENERAL << endl;
  cout << "EMDL_POSTPROCESS_FSC_TRUE = " << EMDL_POSTPROCESS_FSC_TRUE << endl;
  cout << "EMDL_POSTPROCESS_FSC_MASKED = " << EMDL_POSTPROCESS_FSC_MASKED << endl;
  cout << "EMDL_POSTPROCESS_FSC_UNMASKED = " << EMDL_POSTPROCESS_FSC_UNMASKED << endl;
  cout << "EMDL_POSTPROCESS_FSC_RANDOM_MASKED = " << EMDL_POSTPROCESS_FSC_RANDOM_MASKED << endl;
  cout << "EMDL_POSTPROCESS_AMPLCORR_MASKED = " << EMDL_POSTPROCESS_AMPLCORR_MASKED << endl;
  cout << "EMDL_POSTPROCESS_AMPLCORR_UNMASKED = " << EMDL_POSTPROCESS_AMPLCORR_UNMASKED << endl;
  cout << "EMDL_POSTPROCESS_DPR_MASKED = " << EMDL_POSTPROCESS_DPR_MASKED << endl;
  cout << "EMDL_POSTPROCESS_DPR_UNMASKED = " << EMDL_POSTPROCESS_DPR_UNMASKED << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_FIT_CORRELATION = " << EMDL_POSTPROCESS_GUINIER_FIT_CORRELATION << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT = " << EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_FIT_SLOPE = " << EMDL_POSTPROCESS_GUINIER_FIT_SLOPE << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_VALUE_IN = " << EMDL_POSTPROCESS_GUINIER_VALUE_IN << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_VALUE_INVMTF = " << EMDL_POSTPROCESS_GUINIER_VALUE_INVMTF << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_VALUE_WEIGHTED = " << EMDL_POSTPROCESS_GUINIER_VALUE_WEIGHTED << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_VALUE_SHARPENED = " << EMDL_POSTPROCESS_GUINIER_VALUE_SHARPENED << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_VALUE_INTERCEPT = " << EMDL_POSTPROCESS_GUINIER_VALUE_INTERCEPT << endl;
  cout << "EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED = " << EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED << endl;
  cout << "EMDL_POSTPROCESS_MTF_VALUE = " << EMDL_POSTPROCESS_MTF_VALUE << endl;
  cout << "EMDL_POSTPROCESS_RANDOMISE_FROM = " << EMDL_POSTPROCESS_RANDOMISE_FROM << endl;
  cout << "EMDL_POSTPROCESS_UNFIL_HALFMAP1 = " << EMDL_POSTPROCESS_UNFIL_HALFMAP1 << endl;
  cout << "EMDL_POSTPROCESS_UNFIL_HALFMAP2 = " << EMDL_POSTPROCESS_UNFIL_HALFMAP2 << endl;
  cout << "EMDL_SAMPLING_IS_3D = " << EMDL_SAMPLING_IS_3D << endl;
  cout << "EMDL_SAMPLING_IS_3D_TRANS = " << EMDL_SAMPLING_IS_3D_TRANS << endl;
  cout << "EMDL_SAMPLING_HEALPIX_ORDER = " << EMDL_SAMPLING_HEALPIX_ORDER << endl;
  cout << "EMDL_SAMPLING_LIMIT_TILT = " << EMDL_SAMPLING_LIMIT_TILT << endl;
  cout << "EMDL_SAMPLING_OFFSET_RANGE = " << EMDL_SAMPLING_OFFSET_RANGE << endl;
  cout << "EMDL_SAMPLING_OFFSET_STEP = " << EMDL_SAMPLING_OFFSET_STEP << endl;
  cout << "EMDL_SAMPLING_HELICAL_OFFSET_STEP = " << EMDL_SAMPLING_HELICAL_OFFSET_STEP << endl;
  cout << "EMDL_SAMPLING_PERTURB = " << EMDL_SAMPLING_PERTURB << endl;
  cout << "EMDL_SAMPLING_PERTURBATION_FACTOR = " << EMDL_SAMPLING_PERTURBATION_FACTOR << endl;
  cout << "EMDL_SAMPLING_PRIOR_MODE = " << EMDL_SAMPLING_PRIOR_MODE << endl;
  cout << "EMDL_SAMPLING_PSI_STEP = " << EMDL_SAMPLING_PSI_STEP << endl;
  cout << "EMDL_SAMPLING_SIGMA_ROT = " << EMDL_SAMPLING_SIGMA_ROT << endl;
  cout << "EMDL_SAMPLING_SIGMA_TILT = " << EMDL_SAMPLING_SIGMA_TILT << endl;
  cout << "EMDL_SAMPLING_SIGMA_PSI = " << EMDL_SAMPLING_SIGMA_PSI << endl;
  cout << "EMDL_SAMPLING_SYMMETRY = " << EMDL_SAMPLING_SYMMETRY << endl;
  cout << "EMDL_SELECTED = " << EMDL_SELECTED << endl;
  cout << "EMDL_SELECT_PARTICLES_ZSCORE = " << EMDL_SELECT_PARTICLES_ZSCORE << endl;
  cout << "EMDL_SORTED_IDX = " << EMDL_SORTED_IDX << endl;
  cout << "EMDL_STARFILE_MOVIE_PARTICLES = " << EMDL_STARFILE_MOVIE_PARTICLES << endl;
  cout << "EMDL_PERFRAME_CUMULATIVE_WEIGHT = " << EMDL_PERFRAME_CUMULATIVE_WEIGHT << endl;
  cout << "EMDL_PERFRAME_RELATIVE_WEIGHT = " << EMDL_PERFRAME_RELATIVE_WEIGHT << endl;
  cout << "EMDL_RESOLUTION = " << EMDL_RESOLUTION << endl;
  cout << "EMDL_RESOLUTION_ANGSTROM = " << EMDL_RESOLUTION_ANGSTROM << endl;
  cout << "EMDL_RESOLUTION_INVPIXEL = " << EMDL_RESOLUTION_INVPIXEL << endl;
  cout << "EMDL_SPECTRAL_IDX = " << EMDL_SPECTRAL_IDX << endl;
  cout << "EMDL_LAST_LABEL = " << EMDL_LAST_LABEL << endl;

  }*/
