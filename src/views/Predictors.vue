<template lang="pug">
v-layout(style='width: 100%')
  v-navigation-drawer(
    width='420',
    permanent,
    style='position: sticky; height: 100vh',
    absolute=false,
    color='background'
  )
    //- Filters
    .pa-4(
      style='border-radius: var(--v-borderRadius) var(--v-borderRadius) 0 0'
    )
      .text-h5.font-weight-bold(style='color: var(--v-text)') Filters
    v-list.pa-4(
      color='transparent',
      v-for='(filter, filter_i) in filters',
      :key='filter.title + filter_i'
    )
      v-list-item-title.mb-n1.font-weight-bold {{ filter.title }}
      small.v-text.text--darken-2 {{ filter.subtitle }}
      v-list-item-content(v-if='filter.type === "range"')
        div(v-if='data.length')
          v-sparkline(
            :fill='true',
            :gradient='getGradient(filter.min, filter.max, filter.range, filter.step)',
            gradient-direction='left',
            :line-width='2',
            :padding='8',
            :smooth=25,
            :value='filter.tickLabels'
          )
          v-range-slider.mt-n11(
            v-model='filter.range',
            :step='filter.step',
            :max='filter.max',
            :min='filter.min',
            track-color='#d9ebe9',
            track-fill-color='text',
            color='text'
          )
          v-layout.mt-n11
            v-text-field.pa-2(
              v-model='filter.range[0]',
              hide-details,
              filled,
              dense,
              label='min',
              type='number',
              color='text'
            )
            v-text-field.pa-2(
              v-model='filter.range[1]',
              hide-details,
              filled,
              dense,
              label='max',
              type='number',
              color='text'
            )
        v-skeleton-loader.mx-auto(
          v-else,
          type='card, actions',
          max-height='100px'
        )
      v-list-item-content.overflow-x(v-if='filter.type === "chip"')
        div(v-if='data.length')
          v-chip-group(
            v-model='filter.selected',
            mandatory,
          )
            v-tooltip(
              bottom,
              v-for='(item, chip_id) in filter.items',
              :key='item.label + chip_id',
              :disabled='typeof item.description == typeof undefined',
              max-width='275'
            )
              template(v-slot:activator='{ on, attrs }')
                v-chip.pa-4(v-on='on', v-bind='attrs', :outlined='filter.selected == chip_id ? false : true', color='primary') {{ item.label }}
              span {{ item.description }}
        v-layout.text-left(col, v-else)
          v-skeleton-loader.mx-auto(type='chip')
          v-skeleton-loader.mx-auto(type='chip')
      v-list-item-content(v-if='filter.type === "autocomplete"')
        div(v-if='data.length')
          v-autocomplete(
            v-model='filter.selected',
            :items='filter.items',
            auto-select-first,
            chips,
            clearable,
            dense,
            deletable-chips,
            filled,
            multiple,
            :label='filter.title'
          )
        v-layout.text-left(col, v-else)
          v-skeleton-loader.mx-auto(type='card-heading')
      v-divider.mt-4
  v-card.ma-6.ml-md-0(width='100%', height='100%', flat)
    v-card-title.pb-2
      v-col.pb-2
        v-text-field.search(
          v-model='searchVisible',
          dense,
          rounded,
          prepend-inner-icon='mdi-magnify',
          label='Type predictor name, algorithm, author...',
          single-line,
          background-color='accent',
          hide-details,
          filled,
          autofocus,
          color='primary',
          @input='updateSearchReal',
          clearable
        )
        //- Mirror filters in v-chips
        v-expand-transition
          v-card.mx-auto.mt-3.float-left(flat)
            v-chip-group(column)
              TransitionGroup(name='scale-w', mode='in-out', tag='div')
                v-chip(
                  color='primary',
                  v-for='(filter, i) in filterChips',
                  :key='i',
                  clearable,
                  @click='resetFilterByChipId(i)',
                  @click:close='resetFilterByChipId(i)',
                  close
                ) 
                  span.font-weight-light(v-if='filter') {{ filter.title + ': ' }}
                  span.pl-1 {{ getFilterDescription(filter) }}
        v-card.mt-3.float-right(flat, width='200')
          v-select(
            v-model='select',
            clearable,
            outlined,
            dense,
            hide-details,
            label='Sort by',
            :items='headers.filter((item) => item.sortable != false)',
            @input='sortItems'
          )
            template(v-slot:prepend)
              v-icon(@click='flipSortOrder()', v-if='select', color='DarkGray') {{ isSortDescending ? 'mdi-arrow-down' : 'mdi-arrow-up' }}
    v-card-text.pb-2
      v-col
        span(v-if='isDataLoaded && !dataVisible.length') Not found. Try to change search request or filters
        //- id should be id, not name
        PredictorCard.mb-2(
          v-for='(card, card_i) in dataVisible',
          :key='card.predictor + card_i',
          :showSkeleton='card.showSkeleton',
          :title='card.predictor',
          :_id='card.predictor',
          :input='card.input',
          :algorithm_0='card.algorithm_0',
          :algorithm_1='card.algorithm_1',
          :stabilization_sign='card.stabilization_sign',
          :datasets='card.datasets',
          :meta='card.meta',
          :cv='card.cv',
          :mutlitple_point_mutations='card.mutlitple_point_mutations',
          :complexes='card.complexes',
          :hrm_dataset='card.hrm_dataset',
          :hrm_check='card.hrm_check',
          :compared_tools='card.compared_tools',
          :metrics='card.metrics',
          :T='card.T',
          :ph='card.ph',
          :author='card.author',
          :year='card.year',
          :reference='card.reference',
          :download='card.download',
          :server='card.server',
          :doi='card.doi'
        )
</template>

<script lang="ts">
import Vue from 'vue'
import Component from 'vue-class-component'
import PredictorCard from '@/components/PredictorCard.vue'

// import { getPredictors } from '@/utils/api'
import { Predictor } from '@/models/Predictor'
import { getPredictors } from '@/utils/api'
import { namespace } from 'vuex-class'

const AppStore = namespace('AppStore')

@Component({
  props: {},
  components: {
    PredictorCard,
  },
})
export default class Predictors extends Vue {
  @AppStore.Mutation setPredictors!: (predictors: Predictor[]) => void

  searchVisible: String = ''
  searchReal: String | null = null

  select = []

  isSortDescending = true

  isDataLoaded = false

  get dataVisible() {
    // Apply filters
    let result = this.data.filter((item) =>
      this.customFilter('', this.searchReal, item)
    )

    // Values are already sorted in this.sortItems
    // With datasets > 10 000 that approach may affect performance

    // if no data, show skeleton
    if (!this.isDataLoaded && !result.length)
      for (let i = 0; i < 5; i++)
        result.push({
          showSkeleton: true,
          predictor: '',
          input: '',
          algorithm_0: '',
          algorithm_1: '',
          datasets: {
            train: [],
            test: [],
          },
          meta: false,
          cv: '',
          multiple_point_mutations: false,
          complexes: false,
          hrm_dataset: false,
          hrm_check: false,
          compared_tools: '',
          metrics: '',
          T: false,
          ph: true,
          author: '',
          doi: '',
          reference: '',
          download: '',
          server: '',
          year: 2011,
        })
    return result
  }

  updateSearchReal() {
    if (this.searchVisible != '') this.searchReal = this.searchVisible
    else this.searchReal = null
  }

  randn_bm(min, max, skew) {
    let u = 0,
      v = 0
    while (u === 0) u = Math.random() //Converting [0,1) to (0,1)
    while (v === 0) v = Math.random()
    let num = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v)

    num = num / 10.0 + 0.5 // Translate to 0 -> 1
    if (num > 1 || num < 0) num = this.randn_bm(min, max, skew)
    // resample between 0 and 1 if out of range
    else {
      num = Math.pow(num, skew) // Skew
      num *= max - min // Stretch to fill range
      num += min // offset to min
    }
    return num
  }

  tickLabels(start: number, end: number, step: number) {
    // That shuold be pre-computed on backend
    let rand: Array<number> = []
    for (let i = 0; i < step * 1000; i++)
      rand.push(Math.round(this.randn_bm(start, end, 1) * step * 2) / step / 2)
    rand = rand.sort(function (a, b) {
      return a - b
    })
    let randReduced: Map<number, Map<number, number>> = rand.reduce(
      (acc, e) => acc.set(e, (acc.get(e) || 0) + 1),
      new Map()
    )
    return [...randReduced.values()]
  }

  tickLabelsByData(param: string) {
    const values = this.data.map((el) => el[param])
    // Group near values or data-augmentation
    let valuesProcessed
    if (param == 'size') {
      valuesProcessed = values.map((el) => Math.round(el / 50) * 50)
    } else valuesProcessed = values
    const valuesReduced: Map<
      number,
      Map<number, number>
    > = valuesProcessed.reduce(
      (acc, e) => acc.set(e, (acc.get(e) || 0) + 1),
      new Map()
    )
    return [...valuesReduced.values()]
  }

  getGradient(min, max, range, step) {
    let color = this.$vuetify.theme.themes[this.$vuetify.theme.dark ? 'dark' : 'light']['primary']
    let colorDisabled = '#d9ebe9'
    if (isNaN(range[0]) || typeof range[0] === 'string') return [colorDisabled]
    // if (min == -10) console.log(range)

    let gradient: Array<String> = []
    for (let i = min; i < max; i += step) {
      if (i < range[0]) gradient.push(colorDisabled)
      if (i > range[0] && i < range[1]) gradient.push(color as String)
      if (i > range[1]) gradient.push(colorDisabled)
    }

    return gradient
  }

  customFilter(value: any, search: String | null, item: object) {
    // in vue 2.6.9 works only if search string is provided
    // So, we need some magic for user to ignore that bug

    let result = false

    if (search == null) result = true
    else
      result = Object.values(item).some((el) =>
        el.toString().toLowerCase().includes(search.toLowerCase())
      )
    return Object.keys(item).reduce((prev, current) => {
      if (prev === false) return false
      // get suitable filter object
      let currentFilter = this.filters.filter((f) => f.value === current)
      if (currentFilter.length == 0) return true
      // check each type
      if (currentFilter[0].type === 'range')
        if (
          //@ts-ignore
          item[current] < currentFilter[0].range[0] ||
          //@ts-ignore
          item[current] > currentFilter[0].range[1]
        )
          return false
      if (currentFilter[0].type === 'chip') {
        // if we selected 'any'
        if (currentFilter[0].selected == 0) return true
        // process download/server case
        if (
          current == 'download' &&
          'items' in Object.keys(currentFilter[0]) && 
          currentFilter[0].items[currentFilter[0].selected as number]
            .fieldToBe == 'download' &&
          item['download'].length
        )
          return true
        if (
          current == 'download' &&
          'items' in Object.keys(currentFilter[0]) && 
          currentFilter[0].items[currentFilter[0].selected as number]
            .fieldToBe == 'server' &&
          item['server'].length
        )
          return true
        // general
        if ('items' in Object.keys(currentFilter[0]) && 
          item[current] !=
          currentFilter[0].items[currentFilter[0].selected as number].fieldToBe
        )
          return false
      }
      if (currentFilter[0].type === 'autocomplete') {
        // if nothing selected
        if (!(currentFilter[0].selected as string[]).length) return true
        if (
          (currentFilter[0].selected as string[]).some((el) =>
            item[current].includes(el)
          )
        )
          return true
        if (
          !(currentFilter[0].selected as string[]).includes(
            item[current] as string
          )
        )
          return false
      }
      return true
    }, result)
  }

  sortItems(value) {
    let factor = this.isSortDescending ? 1 : -1
    this.data.sort((a, b) => factor * (b[value] - a[value]))
  }

  flipSortOrder() {
    this.isSortDescending = !this.isSortDescending
    if (this.select) this.sortItems(this.select)
  }

  get filterChips() {
    return this.filters.filter((item) => {
      if (item.type === 'range') {
        //@ts-ignore
        if (item.min < item.range[0]) return true
        //@ts-ignore
        if (item.range[1] < item.max) return true
      }
      if (item.type === 'chip') {
        if (item.selected) return true
      }
      if (item.type === 'autocomplete') {
        if ((item.selected as string[]).length) return true
      }
      // TBD for other filter types
      return false
    })
  }

  resetFilterByChipId(chipId) {
    // console.log(this.filterChips[chipId])
    let originalIndex = this.filters.indexOf(this.filterChips[chipId])
    if (this.filters[originalIndex].type === 'range') {
      this.resetRangeSlider(originalIndex)
    }
    if (this.filters[originalIndex].type === 'chip')
      this.filters[originalIndex].selected = 0
    if (this.filters[originalIndex].type === 'autocomplete')
      this.filters[originalIndex].selected = []
  }

  resetRangeSlider(i) {
    // solves reactivity problem
    //@ts-ignore
    Vue.set(this.filters[i].range, 0, this.filters[i].min)
    //@ts-ignore
    Vue.set(this.filters[i].range, 1, this.filters[i].max)
  }

  getFilterDescription(
    filter: typeof this.filters[0] | typeof this.filters[1]
  ) {
    if (filter.type === 'range')
      return (
        //@ts-ignore
        (filter.range[0] == filter.min ? 'min' : filter.range[0]) +
        '→' +
        //@ts-ignore
        (filter.range[1] == filter.max ? 'max' : filter.range[1])
      )
    if (filter.type === 'chip')
      return filter.items[filter.selected as number].label
    //@ts-ignore
    if (filter.type === 'autocomplete') return filter.selected.join(', ')
  }

  created(): void {
    // set correct sidebar size (25 characters)
    this.$vuetify.theme.themes['light'].sidebar_size = '42ch'
    this.$vuetify.theme.themes['dark'].sidebar_size = '42ch'

    document.title = 'Predictors | ' + this.$t('title')
    ;(document.querySelector(
      'meta[name="description"]'
    ) as HTMLElement).setAttribute(
      'content',
      'Predictors | ' + this.$t('title') as string
    )

    getPredictors().then((response) => {
      this.data = response
      this.isDataLoaded = true
      // TODO: setPredictors to AppStore
      this.setPredictors(response)
      const generateAutompleteItems = [
        'algorithm_0',
        'compared_tools',
        'metrics',
        'author',
      ]
      generateAutompleteItems.forEach((fieldName) => {
        // .split(',') and => set => values

        // Creates set of values (authors) or multiple values (of comma separated string)
        if (this.filters.findIndex((el) => el.value == fieldName) == -1) return;
        // @ts-ignore
        this.filters[
          this.filters.findIndex((el) => el.value == fieldName)
        ].items = [
          ...new Set(
            new Set(
              this.data
                .map((el: Predictor) => el[fieldName].split(','))
                .flat(1)
                .map((e) => e.trim())
                .filter((e) => e)
                .sort()
            )
          ),
        ]
      })
      // console.log(this.filters)
    })
  }

  // TODO: hide subtitles inside tooltips
  filters = [
    {
      title: 'Input',
      subtitle: '',
      value: 'input',
      type: 'chip',
      items: [
        { label: 'Any', fieldToBe: undefined, description: undefined },
        { label: 'Sequence', fieldToBe: 'sequence', description: '' },
        { label: 'Structure', fieldToBe: 'structure', description: '' },
      ],
      selected: 0,
    },
    {
      title: 'Algorithm',
      value: 'algorithm_0',
      subtitle: '',
      type: 'autocomplete',
      items: [],
      selected: [],
    },
    {
      title: 'Access',
      subtitle: '',
      value: 'download',
      type: 'chip',
      items: [
        { label: 'Any', fieldToBe: undefined, description: undefined },
        { label: 'Standalone', fieldToBe: 'download', description: '' },
        { label: 'Server', fieldToBe: 'server', description: '' },
      ],
      selected: 0,
    },
    {
      title: 'Meta-predictor',
      subtitle: '',
      value: 'meta',
      type: 'chip',
      items: [
        { label: 'Any', fieldToBe: undefined, description: undefined },
        { label: 'Yes', fieldToBe: true, description: '' },
        { label: 'No', fieldToBe: false, description: '' },
      ],
      selected: 0,
    },
    {
      title: 'Prediction of multiple-point mutations',
      subtitle: '',
      value: 'multiple_point_mutations',
      type: 'chip',
      items: [
        { label: 'Both', fieldToBe: undefined, description: undefined },
        { label: 'Yes', fieldToBe: true, description: '' },
        { label: 'No', fieldToBe: false, description: '' },
      ],
      selected: 0,
    },
    {
      title: 'Setting temperature',
      subtitle: '',
      value: 'T',
      type: 'chip',
      items: [
        { label: 'Any', fieldToBe: undefined, description: undefined },
        { label: 'Yes', fieldToBe: true, description: '' },
        { label: 'No', fieldToBe: false, description: '' },
      ],
      selected: 0,
    },
    {
      title: 'Setting pH',
      subtitle: '',
      value: 'ph',
      type: 'chip',
      items: [
        { label: 'Any', fieldToBe: undefined, description: undefined },
        { label: 'Yes', fieldToBe: true, description: '' },
        { label: 'No', fieldToBe: false, description: '' },
      ],
      selected: 0,
    },
    {
      title: 'Symmetrized training set',
      subtitle: '',
      value: 'hrm_dataset',
      type: 'chip',
      items: [
        { label: 'Any', fieldToBe: undefined, description: undefined },
        { label: 'Yes', fieldToBe: true, description: '' },
        { label: 'No', fieldToBe: false, description: '' },
      ],
      selected: 0,
    },
    {
      title: 'Author',
      value: 'author',
      subtitle: '',
      type: 'autocomplete',
      items: [],
      selected: [],
    },
  ]
  selected = []

  // most of the values is unnecessary :/
  // TODO: rm unused rows
  headers = [
    { text: 'Predictor', align: 'start', sortable: false, value: 'predictor' },
    { text: 'Input', value: 'input' },
    {
      text: 'Algorithm',
      value: 'algorithm_0',
      sortable: false,
      align: 'start',
    },
    { text: '', value: 'algorithm_1', sortable: false, align: 'start' },
    { text: 'Meta', value: 'meta', sortable: false },
    {
      text: 'Multiple-point mutations ',
      value: 'multiple_point_mutations',
      sortable: false,
    },
    {
      text: 'Mutations in protein complexes',
      value: 'complexes',
      sortable: false,
    },
    { text: 'Setting temperature', value: 'T', sortable: false },
    { text: 'Setting pH', value: 'ph', sortable: false },

    { text: 'Year', value: 'year' },
    {
      text: 'Reference',
      value: 'doi',
      sortable: false,
      align: 'start',
      width: '10',
    },
    { text: 'Download', value: 'download', sortable: false, align: 'start' },
    { text: 'Server', value: 'server', sortable: false, align: 'start' },
  ]

  data: Predictor[] = []
}
</script>

<!-- Add "scoped" attribute to limit CSS to this component only -->
<style scoped>
h3 {
  margin: 40px 0 0;
}
ul {
  list-style-type: none;
  padding: 0;
}
li {
  display: inline-block;
  margin: 0 10px;
}

.scale-w-enter,
.scale-w-enter-from,
.scale-w-leave,
.scale-w-leave-to {
  transform: scaleX(0.5);
  /* -webkit-transform: translate(-100px, 0); */
  /* transform: translate(-100px, 0); */
}
.scale-w-enter-active,
.scale-w-leave-active,
.scale-w-move {
  transition: all 0.2s ease;
}
.scale-w-leave-active {
  position: absolute;
}
</style>
