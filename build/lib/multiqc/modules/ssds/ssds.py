#!/usr/bin/env python

""" MultiQC module to parse output from SSDS reports """

from __future__ import print_function
from collections import OrderedDict

import logging
import re
import sys

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

from multiqc.plots import linegraph, bargraph, scatter, table, heatmap, beeswarm

## Import the ssds submodules
#from . import ssstats

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Samtools has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='SSDS QC reports',
        anchor='ssds', target='SSDS QC',
        href='http://genome.cshlp.org/content/early/2012/03/20/gr.130583.111.full.pdf',
        info=" generates statistics for Single Stranded DNA Sequencing (SSDS).")

        # Set up class objects to hold parsed data
        self.sections = list()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

#        # Call submodule functions
#        n['ssstats'] = ssstats.parse_reports(self)
#        if n['ssstats'] > 0:
#            log.info("Found {} stats reports".format(n['ssstats']))

        # Call submodule functions
        n['ssds'] = self.parse_reports()
        if n['ssds'] > 0:
            log.info("Found {} stats reports".format(n['ssds']))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    def parse_reports(self):
        """ Find ssds stats logs and parse their data """

        # Create dictionaries for stats and histogram data
        self.ssds_stats   = dict()
        self.histograms   = dict()
        self.heatMapList  = dict()
        self.frip_headers = OrderedDict()

        self.heatXcats = []
        self.heatYcats = []

        # Loop through all files in folder for SSDS logs
        #for f in sorted(self.find_log_files(config.sp['ssds']['ssstats'])):
        #for f in sorted(self.find_log_files('ssds')):
        for f in self.find_log_files('ssds'):

            # Chop off the file extension (save space on screen)
            reportName = f['s_name']
            reportName = re.sub('\..+','', reportName)

            parsed_data = dict()

            log.info(reportName)
            log.info('Whoop Whoop !!')

            # Loop through the file by line
            for line in f['f'].splitlines():

                #  Split line into columns
                sections = line.split("\t")

                # Lines starting with totinfo are general information
                # This excludes those lines for histogram data (ITR/uH Len ...etc).
                if not line.startswith("totinfo") and not line.startswith("FRIP"):

                    ## Build histograms
                    #  3 column format : type       length   count
                    #  i.e.              ITR_dsDNA  1        100
                    myType =   sections[0]
                    nlen =     sections[1]
                    nCount =   sections[2]

                    # Create a dictionary for this type unless one already exists
                    if not (myType in self.histograms):
                        self.histograms[myType] = dict()

                    # Create a dictionary for this file and type unless one already exists
                    if not (f['fn'] in self.histograms[myType]):
                        self.histograms[myType][f['fn']] = dict()

                    # Add this value : count pair to the histogram
                    self.histograms[myType][f['fn']][int(nlen)] = int(nCount)

                    # Debug message
                    # sys.stderr.write('OK ... ' + myType + ' : ' + nlen + ' : ' + nCount + '\n')

                    continue

                if line.startswith("FRIP"):
                    col1    = sections[0].split("_")
                    reptype = col1[0]
                    dnatype = col1[1]
                    allele  = sections[1]
                    FRIP    = sections[2]

                    if not (dnatype in self.heatMapList):
                        self.heatMapList[dnatype] = dict();

                    if not (reportName in self.heatMapList[dnatype]):
                        self.heatMapList[dnatype][reportName] = dict();

                    if not (allele in self.heatMapList[dnatype][reportName]):
                        self.heatMapList[dnatype][reportName][allele] = dict();

    #               if not (allele in self.heatMapList[f['fn']][species]):
    #                   self.heatMapList[f['fn']][species][allele] = dict();

                    self.frip_headers[allele] = {
                        'title': allele,
                        'description': 'ssDNA (t1) in hotspots / total ssDNA (t1) (%): ' + allele,
                        'max': 100,
                        'min': 0,
                        'suffix': '%',
                        'scale': 'RdYlGn',
                        'format': '{:.2f}%'
                    }

                    if (float(FRIP) < 0.0001):
                        self.heatMapList[dnatype][reportName][allele] = 'NA'
                    else:
                        self.heatMapList[dnatype][reportName][allele] = FRIP

                    self.heatXcats.append(allele)
                    self.heatYcats.append(reportName)

                ## Get Simple Alignment data
                field = sections[1].strip()
                field = field.replace(' ','_')
                value = float(sections[2].strip())
                parsed_data[field] = value

            # If we parsed something
            if len(parsed_data) > 0:

                # Calculate some customized values for our data
                parsed_data['Total_Aligned']   = parsed_data['ssDNA_type1_fragments'] + parsed_data['ssDNA_type2_fragments'] + parsed_data['dsDNA_fragments'] + parsed_data['unclassified_fragments']
                parsed_data['other']           = parsed_data['total_fragments'] - parsed_data['adapter'] - parsed_data['ssDNA_type1_fragments'] - parsed_data['ssDNA_type2_fragments'] - parsed_data['dsDNA_fragments'] - parsed_data['unclassified_fragments']

                if parsed_data['ssDNA_type1_fragments'] > 0:
                    parsed_data['ssDNA_NR']        = (parsed_data['unique_ssDNA_type1_fragments'] / parsed_data['ssDNA_type1_fragments']) * 100
                else:
                    parsed_data['ssDNA_NR']        = 1

                if parsed_data['total_fragments'] > 0:
                    parsed_data['Aligned_percent'] = (parsed_data['Total_Aligned'] / parsed_data['total_fragments']) * 100
                    parsed_data['adapter_percent'] = (parsed_data['adapter'] / parsed_data['total_fragments']) * 100
                else:
                    parsed_data['Aligned_percent'] = 0
                    parsed_data['adapter_percent'] = 0
                    parsed_data['ssDNA_NR']        = 0
                    parsed_data['other']           = 0

                # Work out some percentages for aligned reads
                if 'total_fragments' in parsed_data:

                    # Scan through all keys
                    for k in list(parsed_data.keys()):

                        # Exclude things that have already been calculated
                        if k != 'total_fragments' and k != 'Total_Aligned' and k != 'Aligned_percent'  and k != 'ssDNA_NR' and k != 'other' and k != 'adapter' and k != 'adapter_percent':
                            # Add a percentage for all T1/T2/ds/unc
                            if  parsed_data['Total_Aligned'] > 0:
                                parsed_data['{}_percent'.format(k)] = ((parsed_data[k]+1) / parsed_data['Total_Aligned']) * 100
                            else:
                                parsed_data['{}_percent'.format(k)] = 0

                # Overwrite duplicate
                if f['s_name'] in self.ssds_stats:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))

                self.add_data_source(f)
                self.ssds_stats[reportName] = parsed_data

        # If we have some data to show, carry on
        if len(self.ssds_stats) > 0:

            # Write parsed report data to a file
            self.write_data_file(self.ssds_stats, 'multiqc_ssds_stats')

            ######################################### General Stats Table
            # This is where we populate the general statistics table
            # Each object is a columnar entry
            self.general_stats_headers['total_fragments'] = {
                'title': 'Tot',
                'description': 'Count of 1st end reads in original BAM (millions)',
                'min': 0,
                'max': 1000,
                'scale': 'PuBu',
                'modify': lambda x: x / 1000000,
            }

            #self.general_stats_headers['adapter'] = {
                #'title': 'Ad (Mn)',
                #'description': 'Aligned Fragments (millions)',
                #'min': 0,
                #'max': 1000,
                #'scale': 'PuBu',
                #'modify': lambda x: x / 1000000,
            #}

            self.general_stats_headers['Total_Aligned'] = {
                'title': 'Aln (Mn)',
                'description': 'Aligned Fragments (millions)',
                'min': 0,
                'max': 1000,
                'scale': 'PuBu',
                'modify': lambda x: x / 1000000,
            }

            self.general_stats_headers['Aligned_percent'] = {
                'title': 'Aln (%)',
                'description': 'Aligned Fragments / Total (%)',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.2f}%'
            }

            self.general_stats_headers['adapter_percent'] = {
                'title': 'Adapt (%)',
                'description': 'Adapter / Total reads (%)',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'YlOrRd',
                'format': '{:.2f}%'
            }

            self.general_stats_headers['ssDNA_NR'] = {
                'title': 'NR ssDNA (%)',
                'description': 'Unique / Total type 1 ssDNA (%)',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.2f}%'
            }

            self.general_stats_headers['ssDNA_type1_fragments_percent'] = {
                'title': 'Type-1 (%)',
                'description': 'Type 1 ssDNA (% of aligned)',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.2f}%'
            }

            self.general_stats_headers['ssDNA_type2_fragments_percent'] = {
                'title': 'Type-2 (%)',
                'description': 'Type 2 ssDNA  (% of aligned)',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.2f}%'
            }

            self.general_stats_headers['dsDNA_fragments_percent'] = {
                'title': 'dsDNA (%)',
                'description': 'dsDNA (% of aligned)',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.2f}%'
            }

            self.general_stats_headers['unclassified_fragments_percent'] = {
                'title': 'Unc (%)',
                'description': 'Unclassified  (% of aligned)',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.2f}%'
            }

            ## Don't know what this does, but it does no harm
            for s_name in self.ssds_stats:
                if s_name not in self.general_stats_data:
                    self.general_stats_data[s_name] = dict()
                self.general_stats_data[s_name].update( self.ssds_stats[s_name] )

            ######################################### Begin making sub-plots

            reads = {
                'min': 0,
                'modify': lambda x: float(x) / 1000000.0,
                'decimalPlaces': 2,
                'shared_key': 'read_count'
            }

            ##################################################################################
            ## PLOT 1: Total alignment barplot
            # Specify the order and colors of the different types of aligned / unaligned reads
            b2keys = OrderedDict()
            b2keys['ssDNA_type1_fragments']  =  { 'color': '#437bb1', 'name': 'ssDNA (T1)' }
            b2keys['ssDNA_type2_fragments']  =  { 'color': '#f7a35c', 'name': 'ssDNA (T2)' }
            b2keys['dsDNA_fragments']        =  { 'color': '#e63491', 'name': 'dsDNA' }
            b2keys['unclassified_fragments'] = { 'color': '#b1084c', 'name': 'Unc' }
            b2keys['adapter']                = { 'color': '#7fffd4', 'name': 'Adapter' }
            b2keys['other']                  = { 'color': '#696969', 'name': 'Other' }

            # Configure the total alignment plot
            b2config = {
                'id': 'Total_Alignment_Stats_Plot',
                'title': 'Alignment QC',
                'ylab': '# Fragments',
                'cpswitch_counts_label': 'Number of Fragments'
            }

            # Add the Total alignment stats barplot to the page
            self.add_section( plot = bargraph.plot(self.ssds_stats, b2keys, b2config),
                              name = 'Total Alignment Stats (including unaligned)',
                              anchor = 'total-aln-stats',
                              content = '<p>This module parses the output from <code>SSDS stats</code>.</p>')

            ##################################################################################
            ## PLOT 2: SSDS alignment barplot
            # Specify the order and colors of the different types of aligned reads
            keys = OrderedDict()
            keys['ssDNA_type1_fragments'] =  { 'color': '#437bb1', 'name': 'ssDNA (T1)' }
            keys['ssDNA_type2_fragments'] =  { 'color': '#f7a35c', 'name': 'ssDNA (T2)' }
            keys['dsDNA_fragments'] =        { 'color': '#e63491', 'name': 'dsDNA' }
            keys['unclassified_fragments'] = { 'color': '#b1084c', 'name': 'Unc' }

            # Configure the SSDS Alignment barplot
            pconfig = {
                'id': 'SSDS_Alignment_Stats_Plot',
                'title': 'Types of Aligned SSDS fragments',
                'ylab': '# Fragments',
                'cpswitch_counts_label': 'Number of Fragments'
            }

            # Add the SSDS barplot to the page
            self.add_section( plot = bargraph.plot(self.ssds_stats, keys, pconfig),
                              name = 'SSDS Stats',
                              anchor = 'ssds-alignment-barplot',
                              content = '<p>This module parses the output from <code>SSDS stats</code>.</p>')

            ##################################################################################
            ## PLOT 3: FRIP Heatmap
            # Specify the order and colors of the different types of aligned / unaligned reads

            # Configure the FRiP table
            pconfig = {
                'id': 'frip_stats_table',
                'table_title': 'FRIP Statistics',
                'save_file': True,
                'raw_data_fn':'multiqc_frip_stats'
            }

            self.frip_headers_sorted = sorted(self.frip_headers);

            # Add the Total alignment stats barplot to the page
            for dna in ['ssType1','ssType2','dsDNA','unclassified']:
                sys.stderr.write(dna)

                self.add_section( plot = table.plot(self.heatMapList[dna], self.frip_headers, pconfig),
                                  name = 'FRIP per species (' + dna + ')',
                                  anchor = 'frip-stats-' + dna,
                                  content = '<p>This module parses the output from <code>SSDS stats</code>.</p>')

            ##################################################################################
            ## PLOTs 4-20: SSDS length histograms

            # Loop through all 16 combinations of parameter and type
            # (t1,t2,ds,un) X (frag,ITR,uH,Offset)
            # hard coding the list is the easiest option here (no need to be clever)
            #for hType in ['Fragments_ssDNA_type1','ITR_ssDNA_type1','uH_ssDNA_type1','Offset_ssDNA_type1','Fragments_ssDNA_type2','ITR_ssDNA_type2','uH_ssDNA_type2','Offset_ssDNA_type2','Fragments_dsDNA','ITR_dsDNA','uH_dsDNA','Offset_dsDNA','Fragments_unclassified','ITR_unclassified','uH_unclassified','Offset_unclassified']:
            for hType in ['Fragments_ssDNA_type1','ITR_ssDNA_type1','uH_ssDNA_type1','Offset_ssDNA_type1','Fragments_ssDNA_type2','ITR_ssDNA_type2','uH_ssDNA_type2','Offset_ssDNA_type2','Fragments_dsDNA','ITR_dsDNA','uH_dsDNA','Offset_dsDNA']:

                # If histogram dictionary exists for this combo
                if len(self.histograms[hType]) > 0:

                    ## Make a percentage normalised version of the data
                    data_percent = {}
                    # Debug messages:
                    # sys.stderr.write("Value : " + hType + '\n')
                    # sys.stderr.write("Value : %s" %  self.histograms[hType]['data'].items() + '\n')

                    # Loop through key, value pairs for this histogram
                    for s_name, data in self.histograms[hType].items():

                        # initialize total to 0 (not sure why I do that)
                        total = 0

                        # Create percentage dictionary
                        data_percent[s_name] = OrderedDict()

                        # Debug messages:
                        # sys.stderr.write(s_name + '\n')
                        # sys.stderr.write("Value : %s" %  data  + '\n')

                        # Get total for this histogram
                        total = float( sum( data.values() ) )

                        # Calculate percentages for this histogram
                        for k, v in data.items():
                            if v > 0:
                                data_percent[s_name][k] = (v/total)*100
                            else:
                                data_percent[s_name][k] = 0

                    # Configure histogram plot
                    histConfig = {
                        'id': hType +'_size',
                        'title': hType +' Size',
                        'ylab': 'Fragments (%)',
                        'xlab': hType +' length (bp)',
                        'xDecimals': False,
                        'ymin': 0,
                        'data_labels': [
                            {'name': 'Percent', 'ylab': 'Fragments (%)'},
                            {'name': 'Counts',  'ylab': 'Fragments (count)'},
                        ]
                    }

                    # Add histogram to multi-QC page
                    # - percentage histogram is first (default)
                    # - switch order to reverse
                    self.add_section( plot = linegraph.plot([data_percent, self.histograms[hType]], histConfig),
                                      name = 'SSDS ' + hType,
                                      anchor = 'ssds-stats' + hType,
                                      content = '<p>This module parses the output from <code>SSDS stats</code>.</p>')

        # Return the number of logs that were found
        return len(self.ssds_stats)
