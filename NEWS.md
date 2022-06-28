## ALSroads 0.2.0

* Fix: check and stop if the road computed is a MULTILINESTRING. It means there is a bug somewhere.
* New: new parameter `alsrad_default_parameters$conductivity$sigma_min = 0.1`. The conductivity cannot be lower than this value. Every values less than this value are clamped to this value. This way the conductivity never contain 0s which are impassable value. This change allows for example to cross small unclassified bridges [#57](https://github.com/Jean-Romain/ALSroads/issues/57).
* Fix: shields are computed differently to fix some edge cases with MULTILINESTRING [#55](https://github.com/Jean-Romain/ALSroads/issues/55).
* Fix: [#53](https://github.com/Jean-Romain/ALSroads/issues/53)
