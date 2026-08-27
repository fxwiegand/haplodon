# Changelog

## [0.17.0](https://github.com/fxwiegand/haplodon/compare/v0.16.0...v0.17.0) (2026-08-25)


### Miscellaneous Chores

* release 0.17.0 ([e89fb4f](https://github.com/fxwiegand/haplodon/commit/e89fb4f9a5cf537b0382c1de652f6afa00167d27))

## [0.16.0](https://github.com/fxwiegand/predictosaurus/compare/v0.15.0...v0.16.0) (2026-07-13)


### Features

* Use probabilistic or for score aggregation across haplotypes ([#291](https://github.com/fxwiegand/predictosaurus/issues/291)) ([907ab8c](https://github.com/fxwiegand/predictosaurus/commit/907ab8c22c9617bdddf6b200478ea12820f53b53))

## [0.15.0](https://github.com/fxwiegand/predictosaurus/compare/v0.14.2...v0.15.0) (2026-07-06)


### Features

* Handle start lost consequence in scoring ([#289](https://github.com/fxwiegand/predictosaurus/issues/289)) ([6dfb76c](https://github.com/fxwiegand/predictosaurus/commit/6dfb76c7e73f5ffbf6153943c4c451f8e219bbb8))

## [0.14.2](https://github.com/fxwiegand/predictosaurus/compare/v0.14.1...v0.14.2) (2026-07-02)


### Bug Fixes

* Use 1-based position for GeneBe variant annotation ([7817e91](https://github.com/fxwiegand/predictosaurus/commit/7817e917588d862509ed0550e2a3f56e660a543b))

## [0.14.1](https://github.com/fxwiegand/predictosaurus/compare/v0.14.0...v0.14.1) (2026-06-30)


### Bug Fixes

* Report reference haplotype frequency for every sample ([4e4cc93](https://github.com/fxwiegand/predictosaurus/commit/4e4cc93834bddbeb7752c6f7d6029b991625af2b))

## [0.14.0](https://github.com/fxwiegand/predictosaurus/compare/v0.13.2...v0.14.0) (2026-06-29)


### Features

* Add --report-protein flag ([#280](https://github.com/fxwiegand/predictosaurus/issues/280)) ([1ef6093](https://github.com/fxwiegand/predictosaurus/commit/1ef6093925731b75a3f290e41d36719d789fa140))
* Allow authenticated GeneBe requests via --genebe-email ([1b0be53](https://github.com/fxwiegand/predictosaurus/commit/1b0be533d917277f594cd745f084890820cfa699))

## [0.13.2](https://github.com/fxwiegand/predictosaurus/compare/v0.13.1...v0.13.2) (2026-06-22)


### Bug Fixes

* Classify read support with reciprocal Bayes factor ([#275](https://github.com/fxwiegand/predictosaurus/issues/275)) ([3d93672](https://github.com/fxwiegand/predictosaurus/commit/3d93672787377c0d9f999b6ef4a9cbf736fb3e8d))

## [0.13.1](https://github.com/fxwiegand/predictosaurus/compare/v0.13.0...v0.13.1) (2026-06-14)


### Bug Fixes

* Apply full alt allele and genomic-ordered variants in haplotype protein ([#270](https://github.com/fxwiegand/predictosaurus/issues/270)) ([cf05637](https://github.com/fxwiegand/predictosaurus/commit/cf056379e5ef8546f6e60211688649ecdc2446c5))
* Guard reference codon lookup against positions before the phase ([#271](https://github.com/fxwiegand/predictosaurus/issues/271)) ([27ca08c](https://github.com/fxwiegand/predictosaurus/commit/27ca08c4b48215a8e9676c6d0c67884fa0e0d76a))
* Prevent negative frameshift by using rem_euclid ([#272](https://github.com/fxwiegand/predictosaurus/issues/272)) ([a93aab9](https://github.com/fxwiegand/predictosaurus/commit/a93aab9cae44aa16ab1c71a935a9a504aeb2ed58))

## [0.13.0](https://github.com/fxwiegand/predictosaurus/compare/v0.12.8...v0.13.0) (2026-06-09)


### Features

* Allow writing full hgvsg notation including reference nodes ([#268](https://github.com/fxwiegand/predictosaurus/issues/268)) ([125002d](https://github.com/fxwiegand/predictosaurus/commit/125002db9a10a31d1492a708ab4c8ac9fd2e8bbc))

## [0.12.8](https://github.com/fxwiegand/predictosaurus/compare/v0.12.7...v0.12.8) (2026-06-03)


### Bug Fixes

* Skip nodes with star alternative allele ([#266](https://github.com/fxwiegand/predictosaurus/issues/266)) ([7611d92](https://github.com/fxwiegand/predictosaurus/commit/7611d924624d930782ec375f5e7035bf463047fa))

## [0.12.7](https://github.com/fxwiegand/predictosaurus/compare/v0.12.6...v0.12.7) (2026-05-28)


### Bug Fixes

* Deduplicate haplotypes after intron removal ([#263](https://github.com/fxwiegand/predictosaurus/issues/263)) ([538bd87](https://github.com/fxwiegand/predictosaurus/commit/538bd874f222f218f0fd69f044dac7f56d53dd3f))

## [0.12.6](https://github.com/fxwiegand/predictosaurus/compare/v0.12.5...v0.12.6) (2026-05-27)


### Bug Fixes

* Use shared GeneBears client with lock ([#261](https://github.com/fxwiegand/predictosaurus/issues/261)) ([183ff2b](https://github.com/fxwiegand/predictosaurus/commit/183ff2b52af0457e25cd2852a79f39eb2fa76156))

## [0.12.5](https://github.com/fxwiegand/predictosaurus/compare/v0.12.4...v0.12.5) (2026-05-27)


### Bug Fixes

* Add optional GeneBe cache support ([faa9648](https://github.com/fxwiegand/predictosaurus/commit/faa9648a0866f44796f3d9ec9a537317d0f35abe))

## [0.12.4](https://github.com/fxwiegand/predictosaurus/compare/v0.12.3...v0.12.4) (2026-05-26)


### Bug Fixes

* Prune edges without evidence early ([#258](https://github.com/fxwiegand/predictosaurus/issues/258)) ([d51bfbd](https://github.com/fxwiegand/predictosaurus/commit/d51bfbdded9a4a114b4e7ce07bd017a1e2252636))

## [0.12.3](https://github.com/fxwiegand/predictosaurus/compare/v0.12.2...v0.12.3) (2026-05-11)


### Bug Fixes

* Stop combining CDS paths naively and remove beam search ([#253](https://github.com/fxwiegand/predictosaurus/issues/253)) ([f673507](https://github.com/fxwiegand/predictosaurus/commit/f6735073f3d7d1b793d52ff671b4dc6c5db57fb6))

## [0.12.2](https://github.com/fxwiegand/predictosaurus/compare/v0.12.1...v0.12.2) (2026-05-04)


### Performance Improvements

* Pre-calculate possible node pairs in graph for edge generation ([#251](https://github.com/fxwiegand/predictosaurus/issues/251)) ([2c34149](https://github.com/fxwiegand/predictosaurus/commit/2c341498ade769e83597f1402766da0e8589689f))

## [0.12.1](https://github.com/fxwiegand/predictosaurus/compare/v0.12.0...v0.12.1) (2026-05-03)


### Performance Improvements

* Refactor read support handling and skip variants as early as ([#249](https://github.com/fxwiegand/predictosaurus/issues/249)) ([3d4264c](https://github.com/fxwiegand/predictosaurus/commit/3d4264c108715ef8158229a7a17e2775cb75ab85))

## [0.12.0](https://github.com/fxwiegand/predictosaurus/compare/v0.11.3...v0.12.0) (2026-04-30)


### Features

* Add HGVS notation option and store both forms ([#247](https://github.com/fxwiegand/predictosaurus/issues/247)) ([b0d21c8](https://github.com/fxwiegand/predictosaurus/commit/b0d21c8885610d27e97c1940bdbbc912d3a273f7))

## [0.11.3](https://github.com/fxwiegand/predictosaurus/compare/v0.11.2...v0.11.3) (2026-04-27)


### Bug Fixes

* Skip failed annotations ([#244](https://github.com/fxwiegand/predictosaurus/issues/244)) ([ad30411](https://github.com/fxwiegand/predictosaurus/commit/ad304116ff5c154d19dd1d40bbaad5fb633ea8ff))

## [0.11.2](https://github.com/fxwiegand/predictosaurus/compare/v0.11.1...v0.11.2) (2026-04-24)


### Bug Fixes

* Fix default value for --genome-build parameter ([#242](https://github.com/fxwiegand/predictosaurus/issues/242)) ([1f9e18b](https://github.com/fxwiegand/predictosaurus/commit/1f9e18b2d496d095391f722f7a702cef194dc53a))

## [0.11.1](https://github.com/fxwiegand/predictosaurus/compare/v0.11.0...v0.11.1) (2026-04-23)


### Bug Fixes

* Use genebears with native-tls ([#240](https://github.com/fxwiegand/predictosaurus/issues/240)) ([f2d77f0](https://github.com/fxwiegand/predictosaurus/commit/f2d77f0866bc862559dee25b4dc620920e391c49))

## [0.11.0](https://github.com/fxwiegand/predictosaurus/compare/v0.10.5...v0.11.0) (2026-04-23)


### Features

* Fetch, store and render variant annotations ([#237](https://github.com/fxwiegand/predictosaurus/issues/237)) ([e83356d](https://github.com/fxwiegand/predictosaurus/commit/e83356d10b2ba3b6477013fb2038d3f90d39a6b4))


### Bug Fixes

* Avoid false ref support in multi-allelic sites ([#239](https://github.com/fxwiegand/predictosaurus/issues/239)) ([517a9dd](https://github.com/fxwiegand/predictosaurus/commit/517a9ddc46c06de1ecc16a6bab65b556dc137092))

## [0.10.5](https://github.com/fxwiegand/predictosaurus/compare/v0.10.4...v0.10.5) (2026-04-21)


### Bug Fixes

* Consider all possible start nodes for path enumeration ([#234](https://github.com/fxwiegand/predictosaurus/issues/234)) ([d19f0b0](https://github.com/fxwiegand/predictosaurus/commit/d19f0b0aaedb8d04dd2bf7c4099e05f42483f964))

## [0.10.4](https://github.com/fxwiegand/predictosaurus/compare/v0.10.3...v0.10.4) (2026-04-20)


### Bug Fixes

* Fix index for multiple variants at single genomic position ([#229](https://github.com/fxwiegand/predictosaurus/issues/229)) ([843a255](https://github.com/fxwiegand/predictosaurus/commit/843a25517cdcd2d6022c954badac4c37c23c6740))

## [0.10.3](https://github.com/fxwiegand/predictosaurus/compare/v0.10.2...v0.10.3) (2026-04-20)


### Bug Fixes

* Filter haplotypes lacking supported edges ([#227](https://github.com/fxwiegand/predictosaurus/issues/227)) ([4fec1dd](https://github.com/fxwiegand/predictosaurus/commit/4fec1dd519a772b73dea4eebff8557bc977ecd99))

## [0.10.2](https://github.com/fxwiegand/predictosaurus/compare/v0.10.1...v0.10.2) (2026-04-09)


### Bug Fixes

* Fix supporting reads output properly and revert [#222](https://github.com/fxwiegand/predictosaurus/issues/222) ([#224](https://github.com/fxwiegand/predictosaurus/issues/224)) ([920af54](https://github.com/fxwiegand/predictosaurus/commit/920af547e974e60e4c454e459cad21ccef10fa7d))
* Initialize supporting reads for all samples ([#222](https://github.com/fxwiegand/predictosaurus/issues/222)) ([1b96a0a](https://github.com/fxwiegand/predictosaurus/commit/1b96a0a437d27691b7e4f0335e8d432e2d75599d))

## [0.10.1](https://github.com/fxwiegand/predictosaurus/compare/v0.10.0...v0.10.1) (2026-04-07)


### Bug Fixes

* Fix and refactor graph edge creation and path traversal ([#220](https://github.com/fxwiegand/predictosaurus/issues/220)) ([85ebd92](https://github.com/fxwiegand/predictosaurus/commit/85ebd92a23c45af0b3206b236d3768a6069694de))

## [0.10.0](https://github.com/fxwiegand/predictosaurus/compare/v0.10.0...v0.10.0) (2026-03-06)


### Features

* Add --min-vaf parameter to build subcommand ([#141](https://github.com/fxwiegand/predictosaurus/issues/141)) ([5ca19c3](https://github.com/fxwiegand/predictosaurus/commit/5ca19c3a7c9cad5b03b4f68f3c062d25c887af62))
* Add `--min-prop-present` parameter to build subcommand ([c4542cb](https://github.com/fxwiegand/predictosaurus/commit/c4542cb108d81db203ef1fbb0ab2f55b452ec97f))
* Add amino acid distance metrics and Grantham matrix ([adfb9c1](https://github.com/fxwiegand/predictosaurus/commit/adfb9c11c3bc2dda45144c5673f51c3f09ca991c))
* Add distance matrix and calculation for sneath ([b09058c](https://github.com/fxwiegand/predictosaurus/commit/b09058cb43febfa33a272e627de428bfd6de578b))
* Add haplotype quantification metrics and store per sample likelihoods derived from VAFs of variants on haplotype ([#126](https://github.com/fxwiegand/predictosaurus/issues/126)) ([6b296a7](https://github.com/fxwiegand/predictosaurus/commit/6b296a7a82de3803e85b2fa8397696221ce5023b))
* Add reverse paths for handling features on reverse strand ([897bd81](https://github.com/fxwiegand/predictosaurus/commit/897bd81193fd062b7adae954e373660ecb511607))
* Add score storage and output for transcripts ([d0919f2](https://github.com/fxwiegand/predictosaurus/commit/d0919f2397564b90a66051154097957041157fc8))
* Add subcommand for neoantigen prediction via peptides ([#78](https://github.com/fxwiegand/predictosaurus/issues/78)) ([4523e4f](https://github.com/fxwiegand/predictosaurus/commit/4523e4f0026fb616f3ea51bd37bc311efc02818a))
* Add support for features on backward strand ([e6b2141](https://github.com/fxwiegand/predictosaurus/commit/e6b21411d3fdce1b54a36f115e8cb76d099ea030))
* Allow conversion of haplotype path to protein ([3fa2c75](https://github.com/fxwiegand/predictosaurus/commit/3fa2c7598268cf7f197aa8943f25b992fb80780e))
* Calculate 0-based position of variant on reverse strand ([8321619](https://github.com/fxwiegand/predictosaurus/commit/8321619843be637c57aad43d671701a93d85b971))
* Calculate scores for haplotypes ([38be031](https://github.com/fxwiegand/predictosaurus/commit/38be0312f987f5844f87ff925aa371dbfc5fcafb))
* Calculate scores for haplotypes ([12518fa](https://github.com/fxwiegand/predictosaurus/commit/12518fa1bd7a35677436836a1219f47324caacba))
* Calculate stop penalty ([0aa1914](https://github.com/fxwiegand/predictosaurus/commit/0aa19147fb9300892c7403493aed2afe965420c1))
* Enable serialization and deserialization for Graph ([200b4dd](https://github.com/fxwiegand/predictosaurus/commit/200b4dd8cba77efcd3aef66c0151d0bc845c7488))
* Finish basic score implementation and restructure part of code ([2ba3d76](https://github.com/fxwiegand/predictosaurus/commit/2ba3d76259076d7428422e99a2924ce2f68f4636))
* Implement Epstein distance matrix and lookup for amino acids ([513c0db](https://github.com/fxwiegand/predictosaurus/commit/513c0db004e143b9f9ecc7571dee1236b661c1ae))
* Implement Miyata distance metric and add tests ([4002801](https://github.com/fxwiegand/predictosaurus/commit/40028018d6af78e02773ff5e8a12889b4d16c5cf))
* Implement new score method ([#158](https://github.com/fxwiegand/predictosaurus/issues/158)) ([91f17c0](https://github.com/fxwiegand/predictosaurus/commit/91f17c0f14504c730cab4ed8cbb4581cdabdfb46))
* Implement show command ([ba4ea1d](https://github.com/fxwiegand/predictosaurus/commit/ba4ea1d1776fdab0d180a9cb96d5b1c62fca9122))
* Implement show command ([a2f243e](https://github.com/fxwiegand/predictosaurus/commit/a2f243ecca5cafd9f65724e103acece84b4bfac8))
* Implement show command ([3bc8d30](https://github.com/fxwiegand/predictosaurus/commit/3bc8d3080dae4e5106a9b7113e08fa92c514def3))
* Implement show command for HTML format ([4af776b](https://github.com/fxwiegand/predictosaurus/commit/4af776b69a00b41a73ec8d8e874195d61be291cd))
* Implement show command for tsv format ([87be0d9](https://github.com/fxwiegand/predictosaurus/commit/87be0d984711520da761e534b540fca3f7a3b8f6))
* Implement show command for tsv format ([4616bc6](https://github.com/fxwiegand/predictosaurus/commit/4616bc6b70c4df1916f8ff68301fe5f2756ee249))
* Implement show command for tsv format ([f0a8965](https://github.com/fxwiegand/predictosaurus/commit/f0a89655851eaf15628d247d6a5857c825e4f87a))
* Restructure into multiple subcommands ([#13](https://github.com/fxwiegand/predictosaurus/issues/13)) ([4462a23](https://github.com/fxwiegand/predictosaurus/commit/4462a233f037cdccc87ec5acfc40af5f0f69d458))
* Start working on actual score calculation (WIP) ([24a1262](https://github.com/fxwiegand/predictosaurus/commit/24a1262d419169e8fbe283abc62d7e22f3625c02))
* Support insertions longer 3 ([252d1d6](https://github.com/fxwiegand/predictosaurus/commit/252d1d64577ab4046b0edc4a2795ba3e2d970932))
* Test graph serialization using serde_json ([5074b49](https://github.com/fxwiegand/predictosaurus/commit/5074b49001a8228c67a3c645ae73469812823cec))
* Update phase according to haplotype properly considering coding ([3ae5746](https://github.com/fxwiegand/predictosaurus/commit/3ae5746e7d5a9e08748469c19bbdf3544a67237e))
* Use amino acid distance in score calculation ([70edff7](https://github.com/fxwiegand/predictosaurus/commit/70edff790032fdb75d5652703a343914ecfd4d2f))
* Use beam search for large path enumeration ([#186](https://github.com/fxwiegand/predictosaurus/issues/186)) ([5fb08d6](https://github.com/fxwiegand/predictosaurus/commit/5fb08d62f3bb2f7009939291ca6aabfc57b0e3d8))
* Use duckdb to write graphs ([#35](https://github.com/fxwiegand/predictosaurus/issues/35)) ([19955ea](https://github.com/fxwiegand/predictosaurus/commit/19955eadaac9d3c0f9dadb7480e8265dca71a285))
* Write haplotypes with their contained variants and add Node::hgvs_notation ([#204](https://github.com/fxwiegand/predictosaurus/issues/204)) ([548d3d4](https://github.com/fxwiegand/predictosaurus/commit/548d3d481cb0046a7260aea01d9689d77129da9d))
* Write supporting read information per haplotype ([#214](https://github.com/fxwiegand/predictosaurus/issues/214)) ([9b95451](https://github.com/fxwiegand/predictosaurus/commit/9b95451ddfa873097b26751531122f92f5681980))


### Bug Fixes

* Add max CDS length option and skip long CDSs by default ([#160](https://github.com/fxwiegand/predictosaurus/issues/160)) ([73b7ef5](https://github.com/fxwiegand/predictosaurus/commit/73b7ef5f8b6d795b3f4a84ca1e85d34b5615ee8a))
* Add missing edges ([ae738d6](https://github.com/fxwiegand/predictosaurus/commit/ae738d611de73afb0df99ff68ac46bda486281e3))
* Add mutex to synchronize score writing in parallel loop ([#135](https://github.com/fxwiegand/predictosaurus/issues/135)) ([1f44ae2](https://github.com/fxwiegand/predictosaurus/commit/1f44ae226de3afbb71698b907d76751cc2c991be))
* Adjust predictosaurus for non-unique CDS identifiers ([#60](https://github.com/fxwiegand/predictosaurus/issues/60)) ([e9c54ed](https://github.com/fxwiegand/predictosaurus/commit/e9c54ed58b4d33dd978b5d6ea50d8670da441bed))
* Allow  Ns in reverse_complement ([c684ece](https://github.com/fxwiegand/predictosaurus/commit/c684ece0c47cec7ff0a3f7435c73611bf1ee00e6))
* Cap numerical overshoot when calculating prob_present() ([95180ec](https://github.com/fxwiegand/predictosaurus/commit/95180ecabe51fa4e935d09cb2fec1c92f8d796b7))
* Cap numerical overshoot when calculating prob_present() ([43eefbe](https://github.com/fxwiegand/predictosaurus/commit/43eefbed4014892cad7e3bb69defde516499c815))
* Cap numerical overshoot when calculating prob_present() ([7f2e73c](https://github.com/fxwiegand/predictosaurus/commit/7f2e73c85c8653af2f19a592fd1cf1dc839e8f2f))
* Clip variants to CDS boundaries ([#167](https://github.com/fxwiegand/predictosaurus/issues/167)) ([6895f4e](https://github.com/fxwiegand/predictosaurus/commit/6895f4eb0b9547e631998b5cf20ff1e023d208ef))
* Consider splicing for each transcript by joining paths for CDSs ([#72](https://github.com/fxwiegand/predictosaurus/issues/72)) ([5c766c9](https://github.com/fxwiegand/predictosaurus/commit/5c766c9b1c8b0466da46961dad82d4efd5cc5ac2))
* Convert CDS start to 0-based index ([a5927aa](https://github.com/fxwiegand/predictosaurus/commit/a5927aa91cf7476e929d7912d0a17e794fa08b1b))
* Create output dir if non existent ([79b1383](https://github.com/fxwiegand/predictosaurus/commit/79b1383075071d506bbfbeff13e2622a5295dd54))
* Create path table in bd before writing to it ([c335418](https://github.com/fxwiegand/predictosaurus/commit/c3354189f046ce9a24cd78baa5b5349435a790fe))
* Ensure `vaf` scale domain is fixed between 0 and 1 ([01fd353](https://github.com/fxwiegand/predictosaurus/commit/01fd353a2a86064a9c5d939b1ea2c979e36a1d2a))
* Filter unique paths ([6106aba](https://github.com/fxwiegand/predictosaurus/commit/6106abaf93b85f60d05af801e868573f28f5af5b))
* Fix combinatory extension of CDSs ([#195](https://github.com/fxwiegand/predictosaurus/issues/195)) ([9726647](https://github.com/fxwiegand/predictosaurus/commit/9726647e37faaa4beb2b02def3b912de217935bb))
* Fix comment name ([0ca41b4](https://github.com/fxwiegand/predictosaurus/commit/0ca41b4445feb599fc0e662242ffbc06af2f4b92))
* Fix contains_variant to take target string argument ([#137](https://github.com/fxwiegand/predictosaurus/issues/137)) ([4a1a9e3](https://github.com/fxwiegand/predictosaurus/commit/4a1a9e3d02e48489a261e8c38cf0c8cbfa84f7eb))
* Fix false region slicing subtracting bases via phase ([#192](https://github.com/fxwiegand/predictosaurus/issues/192)) ([80f0ebd](https://github.com/fxwiegand/predictosaurus/commit/80f0ebd8e7322716fdc41a74bd3ee517dbfdac8a))
* Fix feature graph generation ([283feec](https://github.com/fxwiegand/predictosaurus/commit/283feec842e031a288dacbdbd1e28fd0a4ebfb5d))
* Fix impact color scheme ([b8bed5e](https://github.com/fxwiegand/predictosaurus/commit/b8bed5e733d49863efead1100ad1bb18c15755b9))
* Fix impact for amino acid altering variants ([38548a0](https://github.com/fxwiegand/predictosaurus/commit/38548a0e35ecfa95db48069a0c36012b61365abb))
* Fix multiple output file creation ([5ec7a40](https://github.com/fxwiegand/predictosaurus/commit/5ec7a40d6faef16745523d309927aa2f374c29e7))
* Fix out dir creation ([19dad8d](https://github.com/fxwiegand/predictosaurus/commit/19dad8dfec95fc7d742f645551a4b4c9860b60a3))
* Fix path generation ([862dd5d](https://github.com/fxwiegand/predictosaurus/commit/862dd5dd453185e561bcf7c964e334dadba58ac6))
* Fix path generation ([f440d26](https://github.com/fxwiegand/predictosaurus/commit/f440d262413ffea7022e2ab7c570842986674207))
* Fix weights concatenation ([31d9a12](https://github.com/fxwiegand/predictosaurus/commit/31d9a12fe03d7fd300295b229f502711027fde34))
* Fix wrong interpretation of probabilities ([7aeadb4](https://github.com/fxwiegand/predictosaurus/commit/7aeadb4a81462c6b1afd8ec8fba5506d39ac8f23))
* Gracefully handle ambiguous nucleotide 'N' for amino acid changes ([#133](https://github.com/fxwiegand/predictosaurus/issues/133)) ([f0e72a8](https://github.com/fxwiegand/predictosaurus/commit/f0e72a8cd5940e9bed4806014f2f0615199d67ff))
* Handle stop codons in amino acid distance metrics ([b0913bf](https://github.com/fxwiegand/predictosaurus/commit/b0913bf41dde2ac2c5983a91af7d18a46664e640))
* Handle variants at start and end of reference ([794a7e0](https://github.com/fxwiegand/predictosaurus/commit/794a7e081eb1775cb1b9c83639c080980efda6b3))
* Ignore samples present in calls file that are not present in given observations ([2bff91a](https://github.com/fxwiegand/predictosaurus/commit/2bff91a88314031c99d811b46e87d429da52702f))
* Improve error messages ([c00a01f](https://github.com/fxwiegand/predictosaurus/commit/c00a01f64d205d53c8837b8381da758c8552199b))
* Improve error messages ([3f41035](https://github.com/fxwiegand/predictosaurus/commit/3f41035a54c39d093eea1c6ee52f2362ec8006df))
* Make position_in_transcript strand-aware ([#209](https://github.com/fxwiegand/predictosaurus/issues/209)) ([6562c3b](https://github.com/fxwiegand/predictosaurus/commit/6562c3b78fc82157d5b788bf04dba009da1105e9))
* Penalize stop codons instead of methionine in scoring ([3cf8c29](https://github.com/fxwiegand/predictosaurus/commit/3cf8c29d197b44908231e95dba5683937ed20bc2))
* Prevent cycles in forward traversal by checking path containment ([#136](https://github.com/fxwiegand/predictosaurus/issues/136)) ([31f0cf4](https://github.com/fxwiegand/predictosaurus/commit/31f0cf4c3e2b1f4e319cd3d1efc3599a6a2a976f))
* Properly treat CDS/GFF end as inclusive ([#211](https://github.com/fxwiegand/predictosaurus/issues/211)) ([b4bed5d](https://github.com/fxwiegand/predictosaurus/commit/b4bed5ddd6133d4e8aa66348a37352c3fe186517))
* Prune haplotype paths when extending them via more CDSs ([#200](https://github.com/fxwiegand/predictosaurus/issues/200)) ([bc7ba8b](https://github.com/fxwiegand/predictosaurus/commit/bc7ba8b71b0316335e5aa882a6acfb32b29d065b))
* Reduce beam search threshold and beam width ([e295c9b](https://github.com/fxwiegand/predictosaurus/commit/e295c9b01103ca786a3d4bbf9b032fac8615c767))
* Refactor field names ([2378127](https://github.com/fxwiegand/predictosaurus/commit/2378127b45be953b96f2f61aef799bf6a8e2f6fa))
* Refactor get_targets to return only chromosomes that are actually ([#132](https://github.com/fxwiegand/predictosaurus/issues/132)) ([1e1c57b](https://github.com/fxwiegand/predictosaurus/commit/1e1c57b976a7e6c21f7708f226d72cefd36f2df6))
* Refine behaviour for invalid EventProbs ([c1a300d](https://github.com/fxwiegand/predictosaurus/commit/c1a300d7ba7f4462e2774a869e31befd0d945c3b))
* Refine behaviour for invalid EventProbs ([4342223](https://github.com/fxwiegand/predictosaurus/commit/4342223da59defd6eb39c76fce4225765dc3efa9))
* Refine behaviour for invalid EventProbs ([ac541b1](https://github.com/fxwiegand/predictosaurus/commit/ac541b1382f5304913c8aae762e663a96e63151d))
* Refine behaviour for invalid EventProbs ([d06b7c4](https://github.com/fxwiegand/predictosaurus/commit/d06b7c45849c604f7dbb9704ddc89290b9c08620))
* Refine behaviour for invalid EventProbs ([89ae172](https://github.com/fxwiegand/predictosaurus/commit/89ae172f7948ff94469d9a0f7dac4b6615690f72))
* Remove duplicate ref nodes for multiple variants with same position ([daf9ee0](https://github.com/fxwiegand/predictosaurus/commit/daf9ee0718050cd4cf9782471f82fa036fd05408))
* Remove edges between nodes with same position ([8879214](https://github.com/fxwiegand/predictosaurus/commit/8879214419cbb4299da7042e0593228e8b38853e))
* Remove edges between nodes with same position. ([1ab8ede](https://github.com/fxwiegand/predictosaurus/commit/1ab8edecb05bfb90f035ab5d6ac6392bc9b9acb8))
* Remove nodes_in_between check for creating edges ([f1535ad](https://github.com/fxwiegand/predictosaurus/commit/f1535ad857f65319d63a47ad975c5dad6095d580))
* Remove primary key constraint on nodes table ([e62ac37](https://github.com/fxwiegand/predictosaurus/commit/e62ac3759bfbc9ae410b6101b82f53485d28be1b))
* Remove unused and potentially breaking placeholder nodes ([1d6f6bf](https://github.com/fxwiegand/predictosaurus/commit/1d6f6bf9377086b180cd4d6cf5f22789c7b6d468))
* Remove unused import of std::ptr::write ([8307d35](https://github.com/fxwiegand/predictosaurus/commit/8307d35b227af363dbb030be44b440e5ef540144))
* Reverse reference sequence as expected by amino acid translations ([590c415](https://github.com/fxwiegand/predictosaurus/commit/590c4159498b6276f801af2e42159601faee2bde))
* Skip path enumeration for large graphs ([04b7c6c](https://github.com/fxwiegand/predictosaurus/commit/04b7c6caf27561f357dfa8254160549296f6b470))
* Use json5 to serialize and deserialize EventProbs ([3330305](https://github.com/fxwiegand/predictosaurus/commit/333030533326cb24acaec31ff776014a60b7cbc9))
* Use proper hgvsc notation for haplotypes ([#206](https://github.com/fxwiegand/predictosaurus/issues/206)) ([c056c01](https://github.com/fxwiegand/predictosaurus/commit/c056c01512f4484f3a088e06c3fc9688c257c26b))
* Write score output to single TSV file instead of per transcript ([#124](https://github.com/fxwiegand/predictosaurus/issues/124)) ([df0773d](https://github.com/fxwiegand/predictosaurus/commit/df0773d3dfe8748754a4530deb2648f760216c08))


### Performance Improvements

* Add indexes for write_graphs tables ([731ea91](https://github.com/fxwiegand/predictosaurus/commit/731ea91120e6ea8aee835aa0a01e02207c566cbf))
* Add threads option and parallelize build and process command ([#130](https://github.com/fxwiegand/predictosaurus/issues/130)) ([8c83f96](https://github.com/fxwiegand/predictosaurus/commit/8c83f96c88e87616755ec46934a044c1db9dcd48))
* Enhance transcript processing with variant detection ([#80](https://github.com/fxwiegand/predictosaurus/issues/80)) ([413db7f](https://github.com/fxwiegand/predictosaurus/commit/413db7f8697b161da019cdfdeae7f6b80525952f))
* Use transaction for writing graphs via DuckDB ([cdc2873](https://github.com/fxwiegand/predictosaurus/commit/cdc28730bd8d82eed60422b43c2139fce2ed8fa4))


### Miscellaneous Chores

* release 0.1.1 ([8551a11](https://github.com/fxwiegand/predictosaurus/commit/8551a110d543be8b84acdfe0f0e0e2692f17ec15))
* release 0.10.0 ([b055ff6](https://github.com/fxwiegand/predictosaurus/commit/b055ff654e1dc09181708f5fedcec312f7a200ca))
* release 0.2.1 ([054d16a](https://github.com/fxwiegand/predictosaurus/commit/054d16af036fb56ec93bb016ca99cc8cc84126fb))
* release 0.2.8 ([1668725](https://github.com/fxwiegand/predictosaurus/commit/1668725766ee90157cb67c3d4c1a30980d0c66bb))

## [0.10.0](https://github.com/fxwiegand/predictosaurus/compare/v0.9.3...v0.10.0) (2026-03-06)


### Features

* Write supporting read information per haplotype ([#214](https://github.com/fxwiegand/predictosaurus/issues/214)) ([9b95451](https://github.com/fxwiegand/predictosaurus/commit/9b95451ddfa873097b26751531122f92f5681980))

## [0.9.3](https://github.com/fxwiegand/predictosaurus/compare/v0.9.2...v0.9.3) (2026-03-04)


### Bug Fixes

* Properly treat CDS/GFF end as inclusive ([#211](https://github.com/fxwiegand/predictosaurus/issues/211)) ([b4bed5d](https://github.com/fxwiegand/predictosaurus/commit/b4bed5ddd6133d4e8aa66348a37352c3fe186517))

## [0.9.2](https://github.com/fxwiegand/predictosaurus/compare/v0.9.1...v0.9.2) (2026-03-03)


### Bug Fixes

* Make position_in_transcript strand-aware ([#209](https://github.com/fxwiegand/predictosaurus/issues/209)) ([6562c3b](https://github.com/fxwiegand/predictosaurus/commit/6562c3b78fc82157d5b788bf04dba009da1105e9))

## [0.9.1](https://github.com/fxwiegand/predictosaurus/compare/v0.9.0...v0.9.1) (2026-03-02)


### Bug Fixes

* Use proper hgvsc notation for haplotypes ([#206](https://github.com/fxwiegand/predictosaurus/issues/206)) ([c056c01](https://github.com/fxwiegand/predictosaurus/commit/c056c01512f4484f3a088e06c3fc9688c257c26b))

## [0.9.0](https://github.com/fxwiegand/predictosaurus/compare/v0.8.4...v0.9.0) (2026-02-26)


### Features

* Write haplotypes with their contained variants and add Node::hgvs_notation ([#204](https://github.com/fxwiegand/predictosaurus/issues/204)) ([548d3d4](https://github.com/fxwiegand/predictosaurus/commit/548d3d481cb0046a7260aea01d9689d77129da9d))

## [0.8.4](https://github.com/fxwiegand/predictosaurus/compare/v0.8.3...v0.8.4) (2026-02-23)


### Bug Fixes

* Prune haplotype paths when extending them via more CDSs ([#200](https://github.com/fxwiegand/predictosaurus/issues/200)) ([bc7ba8b](https://github.com/fxwiegand/predictosaurus/commit/bc7ba8b71b0316335e5aa882a6acfb32b29d065b))

## [0.8.3](https://github.com/fxwiegand/predictosaurus/compare/v0.8.2...v0.8.3) (2026-02-12)


### Bug Fixes

* Reduce beam search threshold and beam width ([e295c9b](https://github.com/fxwiegand/predictosaurus/commit/e295c9b01103ca786a3d4bbf9b032fac8615c767))

## [0.8.2](https://github.com/fxwiegand/predictosaurus/compare/v0.8.1...v0.8.2) (2026-02-11)


### Bug Fixes

* Fix combinatory extension of CDSs ([#195](https://github.com/fxwiegand/predictosaurus/issues/195)) ([9726647](https://github.com/fxwiegand/predictosaurus/commit/9726647e37faaa4beb2b02def3b912de217935bb))

## [0.8.1](https://github.com/fxwiegand/predictosaurus/compare/v0.8.0...v0.8.1) (2026-02-10)


### Bug Fixes

* Fix false region slicing subtracting bases via phase ([#192](https://github.com/fxwiegand/predictosaurus/issues/192)) ([80f0ebd](https://github.com/fxwiegand/predictosaurus/commit/80f0ebd8e7322716fdc41a74bd3ee517dbfdac8a))

## [0.8.0](https://github.com/fxwiegand/predictosaurus/compare/v0.7.4...v0.8.0) (2026-02-04)


### Features

* Use beam search for large path enumeration ([#186](https://github.com/fxwiegand/predictosaurus/issues/186)) ([5fb08d6](https://github.com/fxwiegand/predictosaurus/commit/5fb08d62f3bb2f7009939291ca6aabfc57b0e3d8))

## [0.7.4](https://github.com/fxwiegand/predictosaurus/compare/v0.7.3...v0.7.4) (2025-12-08)


### Bug Fixes

* Skip path enumeration for large graphs ([04b7c6c](https://github.com/fxwiegand/predictosaurus/commit/04b7c6caf27561f357dfa8254160549296f6b470))

## [0.7.3](https://github.com/fxwiegand/predictosaurus/compare/v0.7.2...v0.7.3) (2025-12-02)


### Bug Fixes

* Clip variants to CDS boundaries ([#167](https://github.com/fxwiegand/predictosaurus/issues/167)) ([6895f4e](https://github.com/fxwiegand/predictosaurus/commit/6895f4eb0b9547e631998b5cf20ff1e023d208ef))

## [0.7.2](https://github.com/fxwiegand/predictosaurus/compare/v0.7.1...v0.7.2) (2025-12-01)


### Bug Fixes

* Convert CDS start to 0-based index ([a5927aa](https://github.com/fxwiegand/predictosaurus/commit/a5927aa91cf7476e929d7912d0a17e794fa08b1b))

## [0.7.1](https://github.com/fxwiegand/predictosaurus/compare/v0.7.0...v0.7.1) (2025-12-01)


### Bug Fixes

* Add max CDS length option and skip long CDSs by default ([#160](https://github.com/fxwiegand/predictosaurus/issues/160)) ([73b7ef5](https://github.com/fxwiegand/predictosaurus/commit/73b7ef5f8b6d795b3f4a84ca1e85d34b5615ee8a))

## [0.7.0](https://github.com/fxwiegand/predictosaurus/compare/v0.6.0...v0.7.0) (2025-11-27)


### Features

* Implement new score method ([#158](https://github.com/fxwiegand/predictosaurus/issues/158)) ([91f17c0](https://github.com/fxwiegand/predictosaurus/commit/91f17c0f14504c730cab4ed8cbb4581cdabdfb46))

## [0.6.0](https://github.com/fxwiegand/predictosaurus/compare/v0.5.0...v0.6.0) (2025-09-07)


### Features

* Add haplotype quantification metrics and store per sample likelihoods derived from VAFs of variants on haplotype ([#126](https://github.com/fxwiegand/predictosaurus/issues/126)) ([6b296a7](https://github.com/fxwiegand/predictosaurus/commit/6b296a7a82de3803e85b2fa8397696221ce5023b))

## [0.5.0](https://github.com/fxwiegand/predictosaurus/compare/v0.4.2...v0.5.0) (2025-09-02)


### Features

* Add --min-vaf parameter to build subcommand ([#141](https://github.com/fxwiegand/predictosaurus/issues/141)) ([5ca19c3](https://github.com/fxwiegand/predictosaurus/commit/5ca19c3a7c9cad5b03b4f68f3c062d25c887af62))


### Bug Fixes

* Add mutex to synchronize score writing in parallel loop ([#135](https://github.com/fxwiegand/predictosaurus/issues/135)) ([1f44ae2](https://github.com/fxwiegand/predictosaurus/commit/1f44ae226de3afbb71698b907d76751cc2c991be))
* Fix contains_variant to take target string argument ([#137](https://github.com/fxwiegand/predictosaurus/issues/137)) ([4a1a9e3](https://github.com/fxwiegand/predictosaurus/commit/4a1a9e3d02e48489a261e8c38cf0c8cbfa84f7eb))
* Gracefully handle ambiguous nucleotide 'N' for amino acid changes ([#133](https://github.com/fxwiegand/predictosaurus/issues/133)) ([f0e72a8](https://github.com/fxwiegand/predictosaurus/commit/f0e72a8cd5940e9bed4806014f2f0615199d67ff))
* Prevent cycles in forward traversal by checking path containment ([#136](https://github.com/fxwiegand/predictosaurus/issues/136)) ([31f0cf4](https://github.com/fxwiegand/predictosaurus/commit/31f0cf4c3e2b1f4e319cd3d1efc3599a6a2a976f))


### Performance Improvements

* Use transaction for writing graphs via DuckDB ([cdc2873](https://github.com/fxwiegand/predictosaurus/commit/cdc28730bd8d82eed60422b43c2139fce2ed8fa4))

## [0.4.2](https://github.com/fxwiegand/predictosaurus/compare/v0.4.1...v0.4.2) (2025-08-20)


### Bug Fixes

* Refactor get_targets to return only chromosomes that are actually ([#132](https://github.com/fxwiegand/predictosaurus/issues/132)) ([1e1c57b](https://github.com/fxwiegand/predictosaurus/commit/1e1c57b976a7e6c21f7708f226d72cefd36f2df6))


### Performance Improvements

* Add threads option and parallelize build and process command ([#130](https://github.com/fxwiegand/predictosaurus/issues/130)) ([8c83f96](https://github.com/fxwiegand/predictosaurus/commit/8c83f96c88e87616755ec46934a044c1db9dcd48))

## [0.4.1](https://github.com/fxwiegand/predictosaurus/compare/v0.4.0...v0.4.1) (2025-08-18)


### Bug Fixes

* Write score output to single TSV file instead of per transcript ([#124](https://github.com/fxwiegand/predictosaurus/issues/124)) ([df0773d](https://github.com/fxwiegand/predictosaurus/commit/df0773d3dfe8748754a4530deb2648f760216c08))

## [0.4.0](https://github.com/fxwiegand/predictosaurus/compare/v0.3.2...v0.4.0) (2025-08-06)


### Features

* Add amino acid distance metrics and Grantham matrix ([adfb9c1](https://github.com/fxwiegand/predictosaurus/commit/adfb9c11c3bc2dda45144c5673f51c3f09ca991c))
* Add distance matrix and calculation for sneath ([b09058c](https://github.com/fxwiegand/predictosaurus/commit/b09058cb43febfa33a272e627de428bfd6de578b))
* Add score storage and output for transcripts ([d0919f2](https://github.com/fxwiegand/predictosaurus/commit/d0919f2397564b90a66051154097957041157fc8))
* Calculate scores for haplotypes ([38be031](https://github.com/fxwiegand/predictosaurus/commit/38be0312f987f5844f87ff925aa371dbfc5fcafb))
* Calculate scores for haplotypes ([12518fa](https://github.com/fxwiegand/predictosaurus/commit/12518fa1bd7a35677436836a1219f47324caacba))
* Calculate stop penalty ([0aa1914](https://github.com/fxwiegand/predictosaurus/commit/0aa19147fb9300892c7403493aed2afe965420c1))
* Finish basic score implementation and restructure part of code ([2ba3d76](https://github.com/fxwiegand/predictosaurus/commit/2ba3d76259076d7428422e99a2924ce2f68f4636))
* Implement Epstein distance matrix and lookup for amino acids ([513c0db](https://github.com/fxwiegand/predictosaurus/commit/513c0db004e143b9f9ecc7571dee1236b661c1ae))
* Implement Miyata distance metric and add tests ([4002801](https://github.com/fxwiegand/predictosaurus/commit/40028018d6af78e02773ff5e8a12889b4d16c5cf))
* Start working on actual score calculation (WIP) ([24a1262](https://github.com/fxwiegand/predictosaurus/commit/24a1262d419169e8fbe283abc62d7e22f3625c02))
* Update phase according to haplotype properly considering coding ([3ae5746](https://github.com/fxwiegand/predictosaurus/commit/3ae5746e7d5a9e08748469c19bbdf3544a67237e))
* Use amino acid distance in score calculation ([70edff7](https://github.com/fxwiegand/predictosaurus/commit/70edff790032fdb75d5652703a343914ecfd4d2f))


### Bug Fixes

* Fix comment name ([0ca41b4](https://github.com/fxwiegand/predictosaurus/commit/0ca41b4445feb599fc0e662242ffbc06af2f4b92))
* Handle stop codons in amino acid distance metrics ([b0913bf](https://github.com/fxwiegand/predictosaurus/commit/b0913bf41dde2ac2c5983a91af7d18a46664e640))
* Penalize stop codons instead of methionine in scoring ([3cf8c29](https://github.com/fxwiegand/predictosaurus/commit/3cf8c29d197b44908231e95dba5683937ed20bc2))
* Remove unused import of std::ptr::write ([8307d35](https://github.com/fxwiegand/predictosaurus/commit/8307d35b227af363dbb030be44b440e5ef540144))

## [0.3.2](https://github.com/fxwiegand/predictosaurus/compare/v0.3.1...v0.3.2) (2025-04-29)


### Bug Fixes

* Remove nodes_in_between check for creating edges ([f1535ad](https://github.com/fxwiegand/predictosaurus/commit/f1535ad857f65319d63a47ad975c5dad6095d580))

## [0.3.1](https://github.com/fxwiegand/predictosaurus/compare/v0.3.0...v0.3.1) (2025-04-29)


### Bug Fixes

* Filter unique paths ([6106aba](https://github.com/fxwiegand/predictosaurus/commit/6106abaf93b85f60d05af801e868573f28f5af5b))

## [0.3.0](https://github.com/fxwiegand/predictosaurus/compare/v0.2.10...v0.3.0) (2025-04-23)


### Features

* Add subcommand for neoantigen prediction via peptides ([#78](https://github.com/fxwiegand/predictosaurus/issues/78)) ([4523e4f](https://github.com/fxwiegand/predictosaurus/commit/4523e4f0026fb616f3ea51bd37bc311efc02818a))


### Bug Fixes

* Cap numerical overshoot when calculating prob_present() ([95180ec](https://github.com/fxwiegand/predictosaurus/commit/95180ecabe51fa4e935d09cb2fec1c92f8d796b7))
* Cap numerical overshoot when calculating prob_present() ([43eefbe](https://github.com/fxwiegand/predictosaurus/commit/43eefbed4014892cad7e3bb69defde516499c815))
* Cap numerical overshoot when calculating prob_present() ([7f2e73c](https://github.com/fxwiegand/predictosaurus/commit/7f2e73c85c8653af2f19a592fd1cf1dc839e8f2f))
* Fix path generation ([862dd5d](https://github.com/fxwiegand/predictosaurus/commit/862dd5dd453185e561bcf7c964e334dadba58ac6))
* Fix path generation ([f440d26](https://github.com/fxwiegand/predictosaurus/commit/f440d262413ffea7022e2ab7c570842986674207))
* Fix weights concatenation ([31d9a12](https://github.com/fxwiegand/predictosaurus/commit/31d9a12fe03d7fd300295b229f502711027fde34))


### Performance Improvements

* Enhance transcript processing with variant detection ([#80](https://github.com/fxwiegand/predictosaurus/issues/80)) ([413db7f](https://github.com/fxwiegand/predictosaurus/commit/413db7f8697b161da019cdfdeae7f6b80525952f))

## [0.2.10](https://github.com/fxwiegand/predictosaurus/compare/v0.2.9...v0.2.10) (2025-02-14)


### Bug Fixes

* Consider splicing for each transcript by joining paths for CDSs ([#72](https://github.com/fxwiegand/predictosaurus/issues/72)) ([5c766c9](https://github.com/fxwiegand/predictosaurus/commit/5c766c9b1c8b0466da46961dad82d4efd5cc5ac2))

## [0.2.9](https://github.com/fxwiegand/predictosaurus/compare/v0.2.8...v0.2.9) (2025-01-17)


### Bug Fixes

* Adjust predictosaurus for non-unique CDS identifiers ([#60](https://github.com/fxwiegand/predictosaurus/issues/60)) ([e9c54ed](https://github.com/fxwiegand/predictosaurus/commit/e9c54ed58b4d33dd978b5d6ea50d8670da441bed))
* Allow  Ns in reverse_complement ([c684ece](https://github.com/fxwiegand/predictosaurus/commit/c684ece0c47cec7ff0a3f7435c73611bf1ee00e6))
* Remove edges between nodes with same position ([8879214](https://github.com/fxwiegand/predictosaurus/commit/8879214419cbb4299da7042e0593228e8b38853e))
* Remove edges between nodes with same position. ([1ab8ede](https://github.com/fxwiegand/predictosaurus/commit/1ab8edecb05bfb90f035ab5d6ac6392bc9b9acb8))
* Remove unused and potentially breaking placeholder nodes ([1d6f6bf](https://github.com/fxwiegand/predictosaurus/commit/1d6f6bf9377086b180cd4d6cf5f22789c7b6d468))
* Reverse reference sequence as expected by amino acid translations ([590c415](https://github.com/fxwiegand/predictosaurus/commit/590c4159498b6276f801af2e42159601faee2bde))

## [0.2.8](https://github.com/fxwiegand/predictosaurus/compare/v0.2.7...v0.2.8) (2025-01-07)


### Miscellaneous Chores

* release 0.2.8 ([1668725](https://github.com/fxwiegand/predictosaurus/commit/1668725766ee90157cb67c3d4c1a30980d0c66bb))

## [0.2.7](https://github.com/fxwiegand/predictosaurus/compare/v0.2.6...v0.2.7) (2025-01-07)


### Bug Fixes

* Ensure `vaf` scale domain is fixed between 0 and 1 ([01fd353](https://github.com/fxwiegand/predictosaurus/commit/01fd353a2a86064a9c5d939b1ea2c979e36a1d2a))

## [0.2.6](https://github.com/fxwiegand/predictosaurus/compare/v0.2.5...v0.2.6) (2024-12-16)


### Bug Fixes

* Create output dir if non existent ([79b1383](https://github.com/fxwiegand/predictosaurus/commit/79b1383075071d506bbfbeff13e2622a5295dd54))
* Fix feature graph generation ([283feec](https://github.com/fxwiegand/predictosaurus/commit/283feec842e031a288dacbdbd1e28fd0a4ebfb5d))

## [0.2.5](https://github.com/fxwiegand/predictosaurus/compare/v0.2.4...v0.2.5) (2024-12-12)


### Bug Fixes

* Remove primary key constraint on nodes table ([e62ac37](https://github.com/fxwiegand/predictosaurus/commit/e62ac3759bfbc9ae410b6101b82f53485d28be1b))

## [0.2.4](https://github.com/fxwiegand/predictosaurus/compare/v0.2.3...v0.2.4) (2024-12-12)


### Bug Fixes

* Fix wrong interpretation of probabilities ([7aeadb4](https://github.com/fxwiegand/predictosaurus/commit/7aeadb4a81462c6b1afd8ec8fba5506d39ac8f23))


### Performance Improvements

* Add indexes for write_graphs tables ([731ea91](https://github.com/fxwiegand/predictosaurus/commit/731ea91120e6ea8aee835aa0a01e02207c566cbf))

## [0.2.3](https://github.com/fxwiegand/predictosaurus/compare/v0.2.2...v0.2.3) (2024-12-10)


### Bug Fixes

* Fix multiple output file creation ([5ec7a40](https://github.com/fxwiegand/predictosaurus/commit/5ec7a40d6faef16745523d309927aa2f374c29e7))

## [0.2.2](https://github.com/fxwiegand/predictosaurus/compare/v0.2.1...v0.2.2) (2024-12-06)


### Bug Fixes

* Improve error messages ([c00a01f](https://github.com/fxwiegand/predictosaurus/commit/c00a01f64d205d53c8837b8381da758c8552199b))
* Improve error messages ([3f41035](https://github.com/fxwiegand/predictosaurus/commit/3f41035a54c39d093eea1c6ee52f2362ec8006df))

## [0.2.1](https://github.com/fxwiegand/predictosaurus/compare/v0.2.0...v0.2.1) (2024-12-06)


### Bug Fixes

* Fix out dir creation ([19dad8d](https://github.com/fxwiegand/predictosaurus/commit/19dad8dfec95fc7d742f645551a4b4c9860b60a3))


### Miscellaneous Chores

* release 0.2.1 ([054d16a](https://github.com/fxwiegand/predictosaurus/commit/054d16af036fb56ec93bb016ca99cc8cc84126fb))

## [0.2.0](https://github.com/fxwiegand/predictosaurus/compare/v0.1.1...v0.2.0) (2024-12-05)


### Features

* Add `--min-prop-present` parameter to build subcommand ([c4542cb](https://github.com/fxwiegand/predictosaurus/commit/c4542cb108d81db203ef1fbb0ab2f55b452ec97f))


### Bug Fixes

* Ignore samples present in calls file that are not present in given observations ([2bff91a](https://github.com/fxwiegand/predictosaurus/commit/2bff91a88314031c99d811b46e87d429da52702f))

## [0.1.1](https://github.com/fxwiegand/predictosaurus/compare/v0.1.0...v0.1.1) (2024-12-04)


### Miscellaneous Chores

* release 0.1.1 ([8551a11](https://github.com/fxwiegand/predictosaurus/commit/8551a110d543be8b84acdfe0f0e0e2692f17ec15))

## 0.1.0 (2024-12-04)


### Features

* Add reverse paths for handling features on reverse strand ([897bd81](https://github.com/fxwiegand/predictosaurus/commit/897bd81193fd062b7adae954e373660ecb511607))
* Add support for features on backward strand ([e6b2141](https://github.com/fxwiegand/predictosaurus/commit/e6b21411d3fdce1b54a36f115e8cb76d099ea030))
* Allow conversion of haplotype path to protein ([3fa2c75](https://github.com/fxwiegand/predictosaurus/commit/3fa2c7598268cf7f197aa8943f25b992fb80780e))
* Calculate 0-based position of variant on reverse strand ([8321619](https://github.com/fxwiegand/predictosaurus/commit/8321619843be637c57aad43d671701a93d85b971))
* Enable serialization and deserialization for Graph ([200b4dd](https://github.com/fxwiegand/predictosaurus/commit/200b4dd8cba77efcd3aef66c0151d0bc845c7488))
* Implement show command ([ba4ea1d](https://github.com/fxwiegand/predictosaurus/commit/ba4ea1d1776fdab0d180a9cb96d5b1c62fca9122))
* Implement show command ([a2f243e](https://github.com/fxwiegand/predictosaurus/commit/a2f243ecca5cafd9f65724e103acece84b4bfac8))
* Implement show command ([3bc8d30](https://github.com/fxwiegand/predictosaurus/commit/3bc8d3080dae4e5106a9b7113e08fa92c514def3))
* Implement show command for HTML format ([4af776b](https://github.com/fxwiegand/predictosaurus/commit/4af776b69a00b41a73ec8d8e874195d61be291cd))
* Implement show command for tsv format ([87be0d9](https://github.com/fxwiegand/predictosaurus/commit/87be0d984711520da761e534b540fca3f7a3b8f6))
* Implement show command for tsv format ([4616bc6](https://github.com/fxwiegand/predictosaurus/commit/4616bc6b70c4df1916f8ff68301fe5f2756ee249))
* Implement show command for tsv format ([f0a8965](https://github.com/fxwiegand/predictosaurus/commit/f0a89655851eaf15628d247d6a5857c825e4f87a))
* Restructure into multiple subcommands ([#13](https://github.com/fxwiegand/predictosaurus/issues/13)) ([4462a23](https://github.com/fxwiegand/predictosaurus/commit/4462a233f037cdccc87ec5acfc40af5f0f69d458))
* Support insertions longer 3 ([252d1d6](https://github.com/fxwiegand/predictosaurus/commit/252d1d64577ab4046b0edc4a2795ba3e2d970932))
* Test graph serialization using serde_json ([5074b49](https://github.com/fxwiegand/predictosaurus/commit/5074b49001a8228c67a3c645ae73469812823cec))
* Use duckdb to write graphs ([#35](https://github.com/fxwiegand/predictosaurus/issues/35)) ([19955ea](https://github.com/fxwiegand/predictosaurus/commit/19955eadaac9d3c0f9dadb7480e8265dca71a285))


### Bug Fixes

* Create path table in bd before writing to it ([c335418](https://github.com/fxwiegand/predictosaurus/commit/c3354189f046ce9a24cd78baa5b5349435a790fe))
* Fix impact color scheme ([b8bed5e](https://github.com/fxwiegand/predictosaurus/commit/b8bed5e733d49863efead1100ad1bb18c15755b9))
* Fix impact for amino acid altering variants ([38548a0](https://github.com/fxwiegand/predictosaurus/commit/38548a0e35ecfa95db48069a0c36012b61365abb))
* Handle variants at start and end of reference ([794a7e0](https://github.com/fxwiegand/predictosaurus/commit/794a7e081eb1775cb1b9c83639c080980efda6b3))
* Refactor field names ([2378127](https://github.com/fxwiegand/predictosaurus/commit/2378127b45be953b96f2f61aef799bf6a8e2f6fa))
* Refine behaviour for invalid EventProbs ([c1a300d](https://github.com/fxwiegand/predictosaurus/commit/c1a300d7ba7f4462e2774a869e31befd0d945c3b))
* Refine behaviour for invalid EventProbs ([4342223](https://github.com/fxwiegand/predictosaurus/commit/4342223da59defd6eb39c76fce4225765dc3efa9))
* Refine behaviour for invalid EventProbs ([ac541b1](https://github.com/fxwiegand/predictosaurus/commit/ac541b1382f5304913c8aae762e663a96e63151d))
* Refine behaviour for invalid EventProbs ([d06b7c4](https://github.com/fxwiegand/predictosaurus/commit/d06b7c45849c604f7dbb9704ddc89290b9c08620))
* Refine behaviour for invalid EventProbs ([89ae172](https://github.com/fxwiegand/predictosaurus/commit/89ae172f7948ff94469d9a0f7dac4b6615690f72))
* Remove duplicate ref nodes for multiple variants with same position ([daf9ee0](https://github.com/fxwiegand/predictosaurus/commit/daf9ee0718050cd4cf9782471f82fa036fd05408))
* Use json5 to serialize and deserialize EventProbs ([3330305](https://github.com/fxwiegand/predictosaurus/commit/333030533326cb24acaec31ff776014a60b7cbc9))
