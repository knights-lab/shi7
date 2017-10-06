#!/usr/bin/env bash
python -m zipapp shi7 -m 'shi7_learner:main' -o shi7_learner.pyz
python -m zipapp shi7 -m 'shi7:main' -o shi7.pyz