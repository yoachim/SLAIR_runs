import numpy as np
import pandas as pd
import sqlite3 as db


def convert_schema(filename, fileout):


    remap_dict ={'RA':'fieldRA', 'dec':'fieldDec', 'mjd':'observationStartMJD', 'exptime':'visitExposureTime', 
    'filter':'filter', 'rotSkyPos':'rotSkyPos', 'nexp':'numExposures',
    'airmass':'airmass', 'FWHMeff':'seeingFwhmEff', 'FWHM_geometric':'seeingFwhmGeom',
    'skybrightness':'skyBrightness', 'night': 'night', 'slewtime':'slewTime', 'fivesigmadepth':'fiveSigmaDepth',
    'alt':'altitude', 'az':'azimuth', 'clouds':'clouds', 'moonAlt':'moonAlt', 'sunAlt':'sunAlt', 'note':'note', 
    'field_id':'fieldId', 'survey_id':'proposalId', 'block_id':'block_id'}


    conn = db.connect(filename)
    df = pd.read_sql('select * from observations;', conn)
    df = df.rename(index=str, columns=remap_dict)
    # Kludge on the visitTime
    df['visitTime'] = 2.*df['numExposures'].values + df['visitExposureTime'].values
    # Dummy column
    df['slewDistance'] = 0.*df['numExposures'].values

    conn.close()
    conn = db.connect(fileout)
    df.to_sql('SummaryAllProps', conn, index=False, if_exists='replace')


if __name__ == '__main__':

    files = {'rolling/rolling_10yrs.db':'rolling/rolling_10yrs_opsim.db',
             'roll_mix/rolling_mix_10yrs.db': 'roll_mix/rolling_mix_10yrs_opsim.db'}

    for infile in files:
        convert_schema(infile, files[infile])

