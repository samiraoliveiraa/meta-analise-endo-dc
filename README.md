# Da disbiose à dor: meta-análise de perfis microbianos associados à endometriose e à dor crônica com aprendizado de máquina

Este repositório contém o pipeline utilizado na meta-análise de estudos de microbiota (amplicon e shotgun) relacionados à endometriose e dor crônica. Os diretórios estão organizados conforme as etapas principais do fluxo de trabalho: obtenção de dados, processamento bioinformático, análises estatísticas e aprendizado de máquina.


### Arquivos:

- ```01 - dados/``` — Contém os scripts relacionados à obtenção das amostras utilizadas na meta-análise.
    - ```ids/``` - arquivos .txt com os IDs das amostras de cada estudo incluído.
    - ```arq_dowload.sh``` - script responsável pelo download automático dos dados utilizando os IDs listados acima.
- ```02 - processamento/``` - Armazena os resultados brutos e intermediários do processamento utilizando ferramentas de bioinformática.
    - ```amplicon/``` - resultados de análises do sequenciamento metataxonômico.
        - ```endo_ata2019/```, ```endo_wei2023/```, ... - diretórios individuais para cada estudo processado.
        - ```meta-analise-int/```, ... - diretórios das análises agrupando os estudos.
        - ```silva-data/``` - contém os arquivos necessários para o treinamento dos classificadores taxonômicos.
    - ```shotgun/``` - resultados de análises do sequenciamento metagenômico.
        - ```etapas/``` - contém os scripts ordenados para o processamento de todos as amostra para os estudos.
- ```03 - análises estatísticas/``` - Inclui notebooks e funções utilizadas nas análises estatísticas e de diversidade microbiana.
    - ```amplicon-plots.ipynb``` - notebook principal de geração de gráficos e análises exploratórias.
    - ```funcoes_est.py``` - arquivo com funções para análises estatísticas.
- ```04 - aprendizado de máquina/``` - Contém scripts e notebooks voltados à identificação de biomarcadores microbianos por meio de algoritmos de *machine learning*.
  - ```batch-effect-correction.r``` - script contendo a correção de *batch effects*.
  - ```taxonomia_amplicon_corrigida.csv``` e ```taxonomia_shotgun_corrigida.tsv``` - tabelas de abundância após a correção de *batch effetcs*.
  - ```modelos_amplicon.py``` e ```modelos_shotgun.py``` - scripts de treinamento e otimização dos modelos.
  - ```modelos-amplicon/``` e ```modelos-shotgun/``` - pastas contendo os modelos treinados.
  - ```biomarcadores.ipynb``` - notebook contendo os resultados e plots dos modelos de aprendizado de máquina.  
