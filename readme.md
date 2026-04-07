# ASM-Lite nf-core Pipeline Documentation

## Pipeline Location

```
/home/zerg/git/asm-lite/nf-core-asmlite/
```

## Structure

```
nf-core-asmlite/
├── main.nf                          # Entry point (routes to bacteria/phages)
├── nextflow.config                  # Parameters, profiles, manifest
├── conf/
│   ├── base.config                  # Resource labels (cpu/memory/time)
│   └── modules.config               # ext.args, ext.prefix, publishDir per process
├── workflows/
│   ├── bacteria.nf                  # Bacteria workflow
│   └── phages.nf                    # Phages workflow
├── subworkflows/local/
│   ├── input_check.nf               # Samplesheet parsing
│   ├── prepare_assembly.nf          # Unicycler -> map reads to assembly
│   ├── mapping.nf                   # Map reads to reference genomes
│   ├── call_variants.nf             # bcftools variant calling
│   └── stats.nf                     # Assembly QC statistics
├── modules/
│   ├── nf-core/                     # Installed nf-core modules (28 modules)
│   └── local/                       # Custom modules
│       ├── contig_abundance.nf
│       └── bracken_add_unclassified.nf
├── bin/
│   └── add_unclassified.sh          # Helper script for Bracken
├── assets/
│   └── multiqc_config.yml
├── tests/
│   ├── samplesheet.csv              # Test samplesheet with per-sample reference FASTA
│   └── references.csv               # Legacy fixture, no longer used for per-sample mapping
└── modules.json                     # nf-core module versions tracking
```

## Usage

```bash
# Bacteria or other non-phage assemblies (Prokka annotation by default)
nextflow run main.nf -profile docker \
  --input samplesheet.csv \
  --kraken2_db /path/to/kraken2_db \
  --outdir results

# Phages (switch annotation to Pharokka)
nextflow run main.nf -profile docker \
  --input samplesheet.csv \
  --pharokka \
  --pharokka_db /path/to/pharokka_db \
  --kraken2_db /path/to/kraken2_db \
  --outdir results

# Same runs without containers, using module `conda` environments
nextflow run main.nf -profile conda \
  --input samplesheet.csv \
  --kraken2_db /path/to/kraken2_db \
  --outdir results

# Optional: override the conda env cache location
export NXF_CONDA_CACHEDIR=/path/to/conda-cache

# Skip flags
--skip_kraken    # Skip Kraken2/Bracken
--skip_variants  # Skip reference mapping + variant calling
--skip_qc        # Skip FastQC
```

### Samplesheet format (`--input`)

```csv
sample,fastq_1,fastq_2,reference
Sample1,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,/path/to/reference1.fna
Sample2,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,/path/to/reference2.fna
```

`reference` is required when variant calling is enabled. To skip reference mapping and VCF calling, use `--skip_variants`.
The old `--references` parameter is now deprecated and ignored if passed.

### Example with your local files

```bash
nextflow run main.nf -profile conda \
  --input /home/zerg/git/nf-core-asmlite/run.csv \
  --kraken2_db /path/to/kraken2_db \
  --outdir results
```

### Runtime profiles

- `-profile docker` for Docker
- `-profile singularity` for Singularity/Apptainer
- `-profile conda` for local Conda environments from module definitions
- `-profile test,<runtime>` for the bundled tiny test dataset, for example `-profile test,conda`

---

## Отличия нового пайплайна от старого

### Контейнеризация

| | Старый | Новый |
|---|--------|-------|
| Runtime | Singularity (.sif файлы) | Docker (BioContainers / Wave) |
| Источник | Локальные сборки (`containers/*.sif`) | Автоматические из nf-core модулей |
| Управление | Ручное | Декларативное в каждом модуле |

### Удалённые шаги

| Шаг | Причина |
|-----|---------|
| **Pilon polishing** | Unicycler в `--mode bold` уже включает полировку. Отдельный Pilon избыточен |
| **SPAdes** | Заменён на Unicycler (который использует SPAdes внутренне) |
| **SnpEff аннотация** | Убрана по решению — можно добавить позже |
| **repair_reads** (BBTools) | Не требуется с fastp |

### Новые возможности

| Возможность | Старый | Новый |
|-------------|--------|-------|
| Phage: Kraken2/Bracken | Нет | Да |
| Phage: Variant calling | Нет | Да |
| CONTIG_ABUNDANCE | Только bacteria | Оба workflow |
| References | Захардкожены в config (r1, r2, r3, r4) | Указываются в `--input` через колонку `reference` |
| Skip flags | Нет | `--skip_kraken`, `--skip_variants`, `--skip_qc` |
| MultiQC | Минимальный | Собирает fastp, fastqc, picard, samtools, quast, bcftools |

### Архитектурные изменения

| | Старый | Новый |
|---|--------|-------|
| Модули | Все кастомные | 28 nf-core + 2 кастомных |
| Версионирование | Нет | Topic-based versioning (nf-core) |
| Config | Всё в nextflow.config | Разделён: base.config + modules.config |
| meta map | `tuple(sample, reads)` | `tuple([id:, single_end:], reads)` — nf-core стандарт |
| Процессы picard | Один `bwa_picard.nf` | Раздельные: BWA_MEM → PICARD_MARKDUPLICATES |
| Вариантный вызов | bcftools mpileup \| call \| view в одном скрипте | Раздельные nf-core модули с промежуточным нормированием |

### Идентичные шаги (логика сохранена)

- FASTP тримминг (`--detect_adapter_for_pe --trim_poly_g`)
- Unicycler сборка (`--mode bold`)
- BWA MEM маппинг (`-M` + Read Groups)
- Kraken2 (`--confidence 0.05 --use-names`) + Bracken (`-r 50 -l S -t 1`)
- bcftools mpileup (`-q 0 -Q 0 --max-depth 10000 -a AD,DP`) + call (`-m --ploidy 1 -v`)
- Фильтрация VCF (`QUAL>=20 && FMT/DP>=10`)
- QUAST, samtools stats/flagstat, Picard WGS/alignment metrics

---

## Как добавлять новые процессы

### Вариант 1: Установить nf-core модуль

nf-core предоставляет 1400+ готовых модулей. Поиск и установка:

```bash
# Активировать окружение nf-core
conda activate nf-core

# Поиск модулей
nf-core modules list remote | grep samtools
nf-core modules info samtools/depth

# Установка
nf-core modules install samtools/depth
```

Модуль появится в `modules/nf-core/samtools/depth/main.nf` и будет автоматически
зарегистрирован в `modules.json`.

**Подключение в subworkflow или workflow:**

```groovy
// 1. Импортировать (с алиасом если нужно)
include { SAMTOOLS_DEPTH } from '../../modules/nf-core/samtools/depth/main'

// 2. Вызвать с правильной сигнатурой (смотри input: в main.nf модуля)
SAMTOOLS_DEPTH(ch_bam_bai, ch_fasta_fai)

// 3. Использовать выходы (смотри output: в main.nf модуля)
ch_depth = SAMTOOLS_DEPTH.out.tsv
```

**Настройка параметров** в `conf/modules.config`:

```groovy
withName: 'SAMTOOLS_DEPTH' {
    ext.args = '-a -q 20'   // аргументы командной строки
    publishDir = [
        path: { "${params.outdir}/depth/${meta.id}" },
        mode: 'copy',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

**Важно:** Перед вызовом модуля всегда читай его `main.nf` — проверь:
- Сколько и каких `input:` ожидает
- Какие `output:` emit доступны
- Какой контейнер используется

### Вариант 2: Создать кастомный модуль

Создай файл в `modules/local/`:

```groovy
// modules/local/my_tool.nf

process MY_TOOL {
    tag "$meta.id"
    label 'process_low'   // ресурсы из base.config

    // Контейнер — укажи Docker образ
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/my_tool:1.0--h1234' :
        'quay.io/biocontainers/my_tool:1.0--h1234' }"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.tsv"), emit: results
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    my_tool \\
        ${args} \\
        --input ${input_file} \\
        --output ${prefix}.results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        my_tool: \$(my_tool --version 2>&1 | head -1)
    END_VERSIONS
    """
}
```

**Подключение:**

```groovy
include { MY_TOOL } from '../../modules/local/my_tool'

MY_TOOL(ch_input)
ch_results = MY_TOOL.out.results
```

### Вариант 3: Создать subworkflow

Создай файл в `subworkflows/local/`:

```groovy
// subworkflows/local/my_analysis.nf

include { TOOL_A } from '../../modules/nf-core/tool_a/main'
include { TOOL_B } from '../../modules/local/tool_b'

workflow MY_ANALYSIS {
    take:
    input_ch   // tuple(meta, file) — опиши формат

    main:
    TOOL_A(input_ch)
    TOOL_B(TOOL_A.out.results)

    emit:
    results = TOOL_B.out.output   // tuple(meta, file)
}
```

**Подключение в workflow:**

```groovy
// workflows/bacteria.nf
include { MY_ANALYSIS } from '../subworkflows/local/my_analysis'

// В main: секции
MY_ANALYSIS(ch_some_input)
ch_multiqc = ch_multiqc.mix(MY_ANALYSIS.out.results.map { meta, f -> f })
```

### Чеклист при добавлении модуля

1. **Проверь сигнатуру** — прочитай `input:` и `output:` в main.nf модуля
2. **Добавь ext.args** в `conf/modules.config` если нужны параметры
3. **Добавь publishDir** в `conf/modules.config` для сохранения результатов
4. **Добавь в MultiQC** если выход совместим — `ch_multiqc.mix(...)`
5. **Проверь контейнер** — `docker pull <image>` перед запуском
6. **Тестируй** — `nextflow run main.nf ... -preview` для проверки синтаксиса

### Полезные команды nf-core

```bash
nf-core modules list local          # Установленные модули
nf-core modules list remote          # Все доступные модули
nf-core modules info <tool>          # Информация о модуле
nf-core modules install <tool>       # Установить модуль
nf-core modules update <tool>        # Обновить модуль
nf-core modules remove <tool>        # Удалить модуль
```

### Доступные labels для ресурсов (conf/base.config)

| Label | CPUs | Memory |
|-------|------|--------|
| process_single | 1 | 4 GB |
| process_low | 2 | 6 GB |
| process_medium | 6 | 8 GB |
| process_high | 12 | 16 GB |

Все ресурсы масштабируются при retry: `* task.attempt`
