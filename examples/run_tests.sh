#!/bin/bash

echo 'Testing PYLAG input... '
echo '>>>> FTLE <<<<<'
python -m MYCOASTLCS -j pylag_ftle.json -i './pylag/pylag_*.nc' -o pylag_.nc
echo '>>>> concentrations and Residence time <<<<<'
python -m MYCOASTLCS -j pylag_conc.json -i './pylag/pylag_*.nc' -o pylag_.nc

echo '>>>>>>> Testing LAGAR input <<<<<<'
echo '>>>> FTLE <<<<<'
python -m MYCOASTLCS -j lagar_ftle.json -i './lagar/*_t_LAGAR.nc' -o lagar_.nc
echo '>>>> Concentrations and Residence time <<<<<'
python -m MYCOASTLCS -j lagar_conc.json -i './lagar/*_t_LAGAR.nc' -o lagar.nc
