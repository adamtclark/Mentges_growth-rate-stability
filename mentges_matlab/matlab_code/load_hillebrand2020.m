%% Load Hillebrand & Kunze 2020 metadata
%
% Loads the metadata (reference, habitat, organism, etc) for all studies
% included in the Hillebrand & Kunze 2020 data set.
%
% This data is taken from Appendix S1 of the original study:
%
%﻿Hillebrand H, Kunze C. 2020. "Meta-analysis on pulse disturbances reveals
% differences in functional and compositional recovery across ecosystems." 
% Ecology Letters 23:575–585.
%
%
% INPUT:
%
% Excel table with column names:
% 
%     {'x'                } counter
%     {'USI'              } study identifier
%     {'resp'             } response variable
%     {'caseID'           } case ID
%     {'studyID'          } study ID
%     {'ref'              } reference
%     {'system'           } system: terrestrial, marine or freshwater
%     {'habitat'          } habitat (wetland, grassland, ..)
%     {'field'            } field or mesocosm
%     {'lat'              } latitude
%     {'long'             } longtitude
%     {'organism'         } organism classes (microbes, plants, ..)
%     {'differentiation'  } difference between experiments (P, crab, ..)
%     {'duration'         } study duration [days]
%     {'sampl'            } How many times sampled
%     {'dist_group'       } disturbance group (fire, removal, ..)
%     {'dist'             } exact disturbance
%     {'dist_type'        } disturbance type: pulse or multipulse
%     {'dist_cat'         } disturbance category (biol, chem, phys, ..)
%     {'dist_magn'        } disturbance magnitude (high, low, mid, NA)
%     {'units_replicates' } 2-35
%     {'open'             }
%     {'func'             } functional or composition
%     {'resp_cat'         } index (univariate diversity; composition, abundance or biomass)
%     {'resil_N'          } number of data points used to calculate slope in resilience measures
%     {'start_LRR'        } log transformed deviance between treatment and control BEFORE the disturbance
%     {'start_var'        } sampling variance of the above 
%     {'T.resist_LRR'       } log transformed deviance between treatment and control AFTER the disturbance 
%     {'resist_var'       }
%     {'T.recov_LRR'        } log transformed deviance between treatment and control AT THE END OF EXPERIMENT
%     {'recov_var'        }
%     {'intcp_rma_day'    }
%     {'se_intcp_rma_day' }
%     {'T.resil_rma_day'    } resilience per day (non-transformed time)
%     {'se_slp_rma_day'   }
%     {'sd_res_rma_day'   }
%     {'T.temp_stab_rma_day'}
%     {'p_rma_day'        }
%     {'intcp_rma_RD'     }
%     {'se_intcp_rma_RD'  }
%     {'T.resil_rma_RD'     } resilience for the entire experiment, time transformed to [0,1]
%     {'se_slp_rma_RD'    }
%     {'sd_res_rma_RD'    }
%     {'T.temp_stab_rma_RD' }
%     {'p_rma_RD'         }
%     {'temporal'         } log(T.temp_stab_rma_RD)
%     {'variab'           }
%     {'resil_var'        }
%     {'col'              }
%     {'size1'            }
%     {'size2'            }
%     {'size3'            }
%     {'size4'            }
%     {'size5'            }
%     {'size6'            }
%     {'temp_var'         }
%     {'org_grp'          }
%     {'unweighted'       }
%     {'resist_var2'      }
%     {'recov_var2'       }
%     {'start_var2'       }
%     {'traject'          }
%     {'traject2'         }
%
%
% Note that:
% ..._day: per day, non-transformed time
% ..._RD: for the entire experiment, for duration [0,1]
%
%
% OUTPUT:
% - Structure "M" containing meta-data of Hillebrand & Kunze 2020 data set,
% with observations in rows and meta-data in columns
% - Table "Table_S1", corresponding to table in supplement, comprising the
% number of observations per organism class in the original data
%
% Table_S1: corresponding to respective table in supplement, containing
% number of observations per organism class in the Hillebrand and Kunze
% 2020 data set.

function [M, Table_S1] = load_hillebrand2020()
%% Import the data

% Import from table
% M = readtable('/Users/am41xite/Documents/Data/Stability data/Hillebrand2020/dataall.xlsx');
M = readtable('dataall.xlsx');

%% Replace NA by NaN

% make sure marker number is not part of the original data
assert(~strcmp(M, '9999.9999'),...
    'Use another marker, this is included in the data!')

input = M.resil_rma_day;
input(strcmp(input, 'NA')) = {'9999.9999'};
output = cellfun(@str2num,input);
output(output==9999.9999) = NaN;
M.resil_rma_day = output;

input = M.resil_rma_RD;
input(strcmp(input, 'NA')) = {'9999.9999'};
output = cellfun(@str2num,input);
output(output==9999.9999) = NaN;
M.resil_rma_RD = output;

input = M.resist_LRR;
input(strcmp(input, 'NA')) = {'9999.9999'};
output = cellfun(@str2num,input);
output(output==9999.9999) = NaN;
M.resist_LRR = output;

%%% Edit 2022/02/12:
%%% This gave an error message with the new Matlab version I'm using now. I
%%% assume it's not needed, as apparently there are no NAs in this column.
% input = M.recov_LRR;
% input(strcmp(input, 'NA')) = {'9999.9999'};
% output = cellfun(@str2num,input);
% output(output==9999.9999) = NaN;
% M.recov_LRR = output;

input = M.temp_stab_rma_RD;
input(strcmp(input, 'NA')) = {'9999.9999'};
output = cellfun(@str2num,input);
output(output==9999.9999) = NaN;
M.temp_stab_rma_RD = output;

input = M.temp_stab_rma_day;
input(strcmp(input, 'NA')) = {'9999.9999'};
output = cellfun(@str2num,input);
output(output==9999.9999) = NaN;
M.temp_stab_rma_day = output;


%% Finer classification of organism classes
% In the original data set, organisms are classified in wide groups, such
% as "vertebrates". I subclassified these here by going back to the
% original studies.

M.organism_original = M.organism;

% Finer classification of category "vertebrates"
M.organism(strcmp(M.ref, 'Bonin && 2011')) = {'birds'};
M.organism(strcmp(M.ref, 'Donohue&&2003')) = [repmat({'macroinvertebrates'},1,4), {'fish'}];
M.organism(strcmp(M.ref, 'Iwata&&2003')) = [{'periphyton'}, repmat({'macroinvertebrates'},1,7), repmat({'fish'},1,6)];
M.organism(strcmp(M.ref, 'Kroll&&2017')) = {'birds'};
M.organism(strcmp(M.ref, 'Syms&Jones2000')) = {'fish'};
M.organism(strcmp(M.ref, 'Williams&&2001')) = {'birds'};

% Finer classification of category "plants"
M.organism(strcmp(M.ref, 'Amami&&2009')) = {'forbs'};
% http://www.edream.ma:8080/jspui/bitstream/123456789/1826/1/Document%20Vegetation%20recolonisation%20of%20a%20Mediterranean%20temporary%20pool%20in%20Morocco%20following%20small-scale%20experimental%20disturbance.pdf
M.organism(strcmp(M.ref, 'Brewer&&1997')) = {'marsh'}; % salt marsh
% https://www.jstor.org/stable/3546601?seq=1
M.organism(strcmp(M.ref, 'Buonopane&&2005')) = {'mixed'}; % non-woody, shrubs
% https://onlinelibrary.wiley.com/doi/full/10.1111/j.0030-1299.2005.13949.x
M.organism(strcmp(M.ref, 'Burge&&2017')) = {'woody'}; % Salix cinerea
% https://onlinelibrary.wiley.com/doi/full/10.1111/avsc.12320
M.organism(strcmp(M.ref, 'Cole&&2003')) = {'mixed'}; % non-woody, shrubs
% https://link.springer.com/article/10.1007%2Fs00267-003-0046-x
M.organism(strcmp(M.ref, 'Cross&Harte2007')) = {'forbs'}; % forbs
% and grasses
% https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/06-1029
M.organism(strcmp(M.ref, 'Flynn&&1995')) = {'marsh'}; % freshwater
% marsh species
% https://link.springer.com/article/10.1007/BF00328426
M.organism(strcmp(M.ref, 'Hester&Mendelssohn2000')) = {'marsh'}; % salt marsh
% https://www.ncbi.nlm.nih.gov/pubmed/11285728
M.organism(strcmp(M.ref, 'Kotanen&&1997')) = {'grasses'}; % grassland
% https://www.jstor.org/stable/2404912?seq=1
M.organism(strcmp(M.ref, 'Kreyling&&2017')) = {'grasses'}; % grassland
% https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12848
M.organism(strcmp(M.ref, 'Longo&&2013')) = {'mixed'}; % grasses and forbs
% https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.12128
M.organism(strcmp(M.ref, 'Montalvo&&1993')) = {'grasses'}; % grassland
% https://onlinelibrary.wiley.com/doi/abs/10.2307/3236107
M.organism(strcmp(M.ref, 'Ogden&&2005')) = {'mixed'}; % forbs and grasses
% https://www.sciencedirect.com/science/article/pii/S0006320705001692
M.organism(strcmp(M.ref, 'Olsen&Klanderud2010')) = {'mixed'}; % shrubs, herbs and bryophytes
% https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.12292
M.organism(strcmp(M.ref, 'Podadera &&2015')) = {'woody'};
% https://www.springerprofessional.de/influence-of-removal-of-a-non-native-tree-species-mimosa-caesalp/5194794
M.organism(strcmp(M.ref, 'Richardson&&2010')) = {'invertebrates'};
% https://link.springer.com/article/10.1007/s10021-010-9317-6
M.organism(strcmp(M.ref, 'Shiels&&2010')) = {'algae'};
% https://www.sciencedirect.com/science/article/pii/S0022098111003133
M.organism(strcmp(M.ref, 'Speed&&2010')) = {'mixed'}; % shrubs, moss, herbs
% https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2745.2010.01685.x
M.organism(strcmp(M.ref, 'Vidra&&2007')) = {'mixed'}; % shrubs, grasses, herbs, vines, ferns
% https://bioone.org/journals/The-Journal-of-the-Torrey-Botanical-Society/volume-134/issue-3/1095-5674(2007)134[410:EOVRON]2.0.CO;2/Effects-of-vegetation-removal-on-native-understory-recovery-in-an/10.3159/1095-5674(2007)134[410:EOVRON]2.0.CO;2.full
M.organism(strcmp(M.ref, 'Wagg&&2017')) = {'grasses'}; % grassland
% https://www.jstor.org/stable/26602229?seq=1#metadata_info_tab_contents
M.organism(strcmp(M.ref, 'Wardle&Jonsson2014')) = {'mixed'}; % understory shrubs, herbs, mosses
% https://esajournals.onlinelibrary.wiley.com/doi/10.1890/13-1666.1
M.organism(strcmp(M.ref, 'Wonkka&&2016')) = {'woody'}; % shrubs
% https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/15-0066

%% Remove all compositional data
% Keep only functional data (biomass, cover, or abundance)

M(strcmp(M.func, 'comp'), :) = [];

%% Clean data
% Remove all organisms with less than 10 observations and all "mixed
% plants"

% Number of observations per organism group
figure('color', 'white', 'position', [680,706,726,272])
hist = histogram(categorical(M.organism));
ylabel('No. observations')

% print rarest organisms
[~,ind] = sort(hist.Values, 'descend');
occurrence = table(hist.Values(ind)', hist.Categories(ind)');
occurrence.Properties.VariableNames = {'num_observations', 'organism_group'};
Table_S1 = occurrence;

% enough data means at least 10 observations
M.enoughData = NaN(size(M.organism));
M.enoughData = ismember(M.organism, hist.Categories(hist.Values>10));

% restrict to organisms with enough data
M(~M.enoughData,:) = [];

% remove the "mixed" plants category (contains both woody and non-woody
% plants that differ strongly in growth rate)
M(strcmp(M.organism, 'mixed'),:) = [];

end




