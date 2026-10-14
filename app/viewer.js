/**
 * Korea forest attribute maps — Earth Engine viewer.
 *
 * Two linked maps under a wipe, so the two product families, or a map and its own uncertainty
 * layer, can be put side by side at the same place and zoom. Clicking reads
 * every published layer at that cell in one request. Controls sit in a panel beside the maps, where
 * the wipe cannot cover them.
 *
 * Paste into the Earth Engine Code Editor and publish as an App.
 *
 * The block below is generated from app/layers.json and data/species_colors.csv. Display ranges are
 * measured per band and pooled across every layer carrying it, so the same colour means the same
 * value on both sides of the wipe. Do not edit it by hand.
 */

// ---- development layout controls ---------------------------------------------------------------
// Keep this set to 'auto' before publishing the Earth Engine App.
var LAYOUT_TEST_MODE = 'auto';
// Allowed values: 'auto', 'mobile', 'desktop'.
var MOBILE_BREAKPOINT = 760;
var DEFAULT_LANGUAGE = 'en';
// Allowed values: 'en', 'ko'. Readers can switch languages without reloading the App.
var CURRENT_LOCATION_ZOOM = 15;
var CURRENT_LOCATION_ICON = 'data:image/svg+xml;base64,' +
  'PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyOCIgaGVpZ2h0PSIyOCIgdmlld0JveD0iMCAwIDI4IDI4' +
  'Ij48cmVjdCB4PSIwLjUiIHk9IjAuNSIgd2lkdGg9IjI3IiBoZWlnaHQ9IjI3IiByeD0iMi41IiBmaWxsPSIjRkZGRkZGIiBzdHJva2U9IiNC' +
  'REJEQkQiLz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSg2IDYpIHNjYWxlKC42NjY2NjY3KSI+PHBhdGggZmlsbD0iIzNDNDA0MyIgZD0iTTEy' +
  'IDhjLTIuMjEgMC00IDEuNzktNCA0czEuNzkgNCA0IDQgNC0xLjc5IDQtNC0xLjc5LTQtNC00em04Ljk0IDNDMjAuNDggNi44MyAxNy4xNyAz' +
  'LjUyIDEzIDMuMDZWMWgtMnYyLjA2QzYuODMgMy41MiAzLjUyIDYuODMgMy4wNiAxMUgxdjJoMi4wNmMuNDYgNC4xNyAzLjc3IDcuNDggNy45' +
  'NCA3Ljk0VjIzaDJ2LTIuMDZjNC4xNy0uNDYgNy40OC0zLjc3IDcuOTQtNy45NEgyM3YtMmgtMi4wNnpNMTIgMTljLTMuODcgMC03LTMuMTMt' +
  'Ny03czMuMTMtNyA3LTcgNyAzLjEzIDcgNy0zLjEzIDctNyA3eiIvPjwvZz48L3N2Zz4=';
var BASEMAP_CONTROL_ICON = 'data:image/svg+xml;base64,' +
  'PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIzMiIgaGVpZ2h0PSIzMiIgdmlld0JveD0iMCAwIDMyIDMy' +
  'Ij48cmVjdCB4PSIwLjUiIHk9IjAuNSIgd2lkdGg9IjMxIiBoZWlnaHQ9IjMxIiByeD0iMi41IiBmaWxsPSIjRkZGRkZGIiBzdHJva2U9IiNC' +
  'REJEQkQiLz48ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSg0IDQpIj48cGF0aCBmaWxsPSIjM0M0MDQzIiBkPSJNMTEuOTkgMTguNTQgNC42MiAx' +
  'Mi44MSAzIDE0LjA3bDkgNyA5LTctMS42My0xLjI3LTcuMzggNS43NHpNMTIgMTZsNy4zNi01LjczTDIxIDlsLTktNy05IDcgMS42MyAxLjI3' +
  'TDEyIDE2em0wLTExLjQ3TDE3Ljc0IDkgMTIgMTMuNDcgNi4yNiA5IDEyIDQuNTN6Ii8+PC9nPjwvc3ZnPg==';
var READOUT_HIGHLIGHT_COLOR = '#000080';
var MOBILE_TITLE_FONT_SIZE = '19px';
var MOBILE_BODY_FONT_SIZE = '16px';
var MOBILE_SMALL_FONT_SIZE = '14px';
var MOBILE_MAP_PORTRAIT_RATIO = 0.76;
var MOBILE_MAP_LANDSCAPE_RATIO = 0.68;

// >>> generated
var VIS = {
  "biomass": {
    "max": 273.3392,
    "min": 60.9964,
    "n_sampled": 9588394,
    "range_from": "p5-p95 pooled over every layer carrying this band",
    "unit": "Mg ha-1"
  },
  "biomass_sd": {
    "max": 82.0048,
    "min": 51.2767,
    "n_sampled": 4794197,
    "range_from": "p5-p95 pooled over every layer carrying this band",
    "unit": "ensemble SD, Mg ha-1"
  },
  "conifer_fraction": {
    "max": 1.0,
    "min": 0.0,
    "range_from": "definitional",
    "unit": "fraction of basal area"
  },
  "conifer_fraction_sd": {
    "max": 0.2903,
    "min": 0.0756,
    "n_sampled": 4794197,
    "range_from": "p5-p95 pooled over every layer carrying this band",
    "unit": "ensemble SD, fraction of basal area"
  },
  "dominant_species_uncertainty": {
    "max": 1.0,
    "min": 0.0,
    "range_from": "definitional",
    "unit": "1 - vote margin"
  },
  "nearest_reference_distance": {
    "max": 0.9217,
    "min": 0.6867,
    "n_sampled": 4794197,
    "range_from": "p5-p95 pooled over every layer carrying this band",
    "unit": "1 - terminal-node proximity"
  },
  "shannon": {
    "max": 1.8211,
    "min": 0.3376,
    "n_sampled": 9588394,
    "range_from": "p5-p95 pooled over every layer carrying this band",
    "unit": "Shannon H\u2032"
  },
  "shannon_sd": {
    "max": 0.5523,
    "min": 0.3351,
    "n_sampled": 4794197,
    "range_from": "p5-p95 pooled over every layer carrying this band",
    "unit": "ensemble SD, Shannon H\u2032"
  }
};
var SPECIES = ["E1E1E1", "E1E1E1", "E1E1E1", "E1E1E1", "C9A227", "FF73DF", "78C679", "E1E1E1", "E1E1E1", "FA9FB5", "267300", "E1E1E1", "E1E1E1", "6A51A3", "E1E1E1", "E1E1E1", "08519C", "E1E1E1", "41B6C4", "9ECAE1", "E1E1E1", "E1E1E1", "1B9E77", "FEC44F", "E1E1E1", "E1E1E1", "A63603", "FDAE6B", "F16913", "E1E1E1", "E1E1E1", "E1E1E1", "E1E1E1", "E1E1E1", "E1E1E1", "E1E1E1", "E1E1E1"];
var SPECIES_MAJOR = [{"code": 4, "hex": "C9A227", "label": "Carpinus laxiflora"}, {"code": 5, "hex": "FF73DF", "label": "Carpinus tschonoskii"}, {"code": 6, "hex": "78C679", "label": "Castanea crenata"}, {"code": 9, "hex": "FA9FB5", "label": "Cornus controversa"}, {"code": 10, "hex": "267300", "label": "Cryptomeria japonica"}, {"code": 13, "hex": "6A51A3", "label": "Larix kaempferi"}, {"code": 16, "hex": "08519C", "label": "Pinus densiflora"}, {"code": 18, "hex": "41B6C4", "label": "Pinus rigida"}, {"code": 19, "hex": "9ECAE1", "label": "Pinus thunbergii"}, {"code": 22, "hex": "1B9E77", "label": "Prunus serrulata var. pubescens"}, {"code": 23, "hex": "FEC44F", "label": "Quercus acutissima"}, {"code": 26, "hex": "A63603", "label": "Quercus mongolica"}, {"code": 27, "hex": "FDAE6B", "label": "Quercus serrata"}, {"code": 28, "hex": "F16913", "label": "Quercus variabilis"}];
var SPECIES_NAME = ["Acer pictum var. mono", "Alnus firma", "Betula pendula", "Betula schmidtii", "Carpinus laxiflora", "Carpinus tschonoskii", "Castanea crenata", "Celtis sinensis", "Chamaecyparis obtusa", "Cornus controversa", "Cryptomeria japonica", "Fraxinus chinensis var. rhynchophylla", "Juglans mandshurica", "Larix kaempferi", "Liriodendron tulipifera", "Morus indica", "Pinus densiflora", "Pinus koraiensis", "Pinus rigida", "Pinus thunbergii", "Platycarya strobilacea", "Populus × tomentiglandulosa", "Prunus serrulata var. pubescens", "Quercus acutissima", "Quercus aliena", "Quercus dentata", "Quercus mongolica", "Quercus serrata", "Quercus variabilis", "Robinia pseudoacacia", "Salix pierotii", "Styrax japonicus", "Tilia amurensis", "Tilia mandshurica", "Ulmus davidiana", "Zelkova serrata", "Other"];
// <<< generated

// Korean common names in published class-code order, from data/species_codes.csv.
var SPECIES_NAME_KO = [
  '고로쇠나무', '사방오리', '자작나무', '박달나무', '서어나무', '개서어나무', '밤나무',
  '팽나무', '편백', '층층나무', '삼나무', '물푸레나무', '가래나무',
  '일본잎갈나무', '백합나무', '산뽕나무', '소나무', '잣나무', '리기다소나무',
  '곰솔', '굴피나무', '은사시나무', '잔털벚나무', '상수리나무', '갈참나무', '떡갈나무',
  '신갈나무', '졸참나무', '굴참나무', '아까시나무', '버드나무', '때죽나무', '피나무',
  '찰피나무', '느릅나무', '느티나무', '기타'
];

// ---- interface language -----------------------------------------------------------------------
// English strings stay as the internal selector keys. Korean is display-only, so changing the
// language cannot change a selected layer, place, map side, centre or zoom.
var UI_TEXT = {
  en: {
    'header.title': 'Korea forest attribute maps',
    'header.author_short': 'Seunguk Kim, SNU',
    'info.author': 'Seunguk Kim, Seoul National University',
    'info.description': 'Four forest attributes over South Korea at 30 m, from the 8th National ' +
      'Forest Inventory (2021–2025), mapped two ways from Sentinel-2, climate, terrain and the ' +
      'national forest-type map.',
    'info.data_link': 'Data (CC BY 4.0)',
    'info.code_link': 'Code (MIT)',
    'info.notes_toggle': 'What to know before using these',
    'map.left_title': 'Left of the wipe',
    'map.right_title': 'Right of the wipe',
    'layout.compare': 'Compare two maps',
    'layout.left_only': 'Left map only',
    'layout.right_only': 'Right map only',
    'settings.family': 'Product family',
    'settings.quantity': 'Quantity',
    'family.direct': 'Direct prediction',
    'family.imputation': 'Plot-index imputation',
    'family.basemap': 'Basemap only',
    'quantity.conifer_fraction': 'Coniferous fraction',
    'quantity.dominant_species': 'Dominant species',
    'quantity.biomass': 'Aboveground biomass',
    'quantity.shannon': 'Tree-species diversity',
    'quantity.conifer_fraction_sd': 'Coniferous fraction — uncertainty',
    'quantity.biomass_sd': 'Aboveground biomass — uncertainty',
    'quantity.shannon_sd': 'Tree-species diversity — uncertainty',
    'quantity.dominant_species_uncertainty': 'Dominant species — uncertainty',
    'quantity.plot_index': 'Matched plot index',
    'quantity.nearest_reference_distance': 'Nearest-reference distance',
    'place.south_korea': 'South Korea',
    'place.jirisan': 'Jirisan National Park',
    'place.jeju': 'Jeju Island',
    'place.lower_geum': 'Lower Geum River basin',
    'mobile.settings_toggle': 'Show / hide map settings',
    'mobile.readout_toggle': 'Show / hide cell readout',
    'mobile.info_toggle': 'Show / hide dataset information',
    'language.toggle': '한국어 / English',
    'basemap.control': 'Basemap settings',
    'basemap.map': 'Map',
    'basemap.satellite': 'Satellite',
    'basemap.labels': 'Show labels',
    'location.button': 'Show my location',
    'location.layer': 'Current location',
    'legend.species_others': 'Others (23 species, 5 % of area)',
    'legend.plot_index': 'colours repeat; they mark where the matched plot changes',
    'range.pooled_p5_p95': 'p5-p95 pooled over every layer carrying this band',
    'range.definitional': 'definitional',
    'unit.biomass': 'Mg ha-1',
    'unit.biomass_sd': 'ensemble SD, Mg ha-1',
    'unit.conifer_fraction': 'fraction of basal area',
    'unit.conifer_fraction_sd': 'ensemble SD, fraction of basal area',
    'unit.dominant_species_uncertainty': '1 - vote margin',
    'unit.nearest_reference_distance': '1 - terminal-node proximity',
    'unit.shannon': 'Shannon H′',
    'unit.shannon_sd': 'ensemble SD, Shannon H′',
    'readout.prompt': 'Click the map to read every layer at that cell.',
    'readout.loading': 'reading…',
    'readout.no_data': 'no data here',
    'readout.outside_extent': 'outside the mapped forest extent',
    'readout.unknown_species': 'Unknown species',
    'readout.other_species': 'Others (pooled)',
    'readout.conifer_fraction': 'Coniferous fraction',
    'readout.biomass': 'Aboveground biomass (Mg ha-1)',
    'readout.shannon': 'Tree-species diversity (Shannon H′)',
    'readout.dominant_species': 'Dominant species',
    'readout.direct': 'Direct',
    'readout.imputed': 'Imputed',
    'readout.matched_plot': 'Matched reference plot',
    'readout.plot_id': 'Plot ID',
    'readout.distance': 'Distance',
    'coordinate.north': 'N',
    'coordinate.east': 'E',
    'note.spatial_pattern': 'Read these as broad spatial patterns, not as field values for a 30 m cell.',
    'note.performance': 'Coniferous fraction and dominant species can be read cell by cell. Biomass ' +
      'and tree-species diversity are much weaker at a single cell and are meant to be read at 1 km ' +
      'or coarser; for those two the plot-index family scores below a mean predictor, so prefer the ' +
      'direct maps there.',
    'note.uncertainty': 'The uncertainty layers rank relative confidence, and the nearest-reference ' +
      'distance says how unlike the cell its matched plot was. Neither is a prediction interval.',
    'note.dominant_species': 'Dominant-species accuracy is concentrated in the common species. ' +
      'Per-class accuracy ships with the dataset in dominant_species_table.csv.',
    'note.plot_index': 'The plot-index family trades per-attribute accuracy for keeping the ' +
      'combination of attributes that occurred together on one real inventory plot. Matching the ' +
      'plot number assigned to each location with the plot number in the NFI-8 dataset allows users ' +
      'to access the full survey record for that plot.'
  },
  ko: {
    'header.title': '한국산림공간정보',
    'header.author_short': '김승욱, 서울대학교',
    'info.author': '김승욱, 서울대학교',
    'info.description': '제8차 국가산림자원조사(2021–2025)와 Sentinel-2·기후·지형·임상도를 바탕으로, ' +
      '2종류의 산출 방식으로 제작한 대한민국 전역 30 m 해상도의 산림 속성 지도입니다.',
    'info.data_link': '데이터 (CC BY 4.0)',
    'info.code_link': '코드 (MIT)',
    'info.notes_toggle': '사용 전 알아둘 점',
    'map.left_title': '왼쪽 지도',
    'map.right_title': '오른쪽 지도',
    'layout.compare': '두 지도 비교',
    'layout.left_only': '왼쪽 지도만',
    'layout.right_only': '오른쪽 지도만',
    'settings.family': '지도 산출 방식',
    'settings.quantity': '표시 항목',
    'family.direct': '직접 예측',
    'family.imputation': '표본점번호 할당',
    'family.basemap': '배경지도만',
    'quantity.conifer_fraction': '침엽수 비율',
    'quantity.dominant_species': '우점 수종',
    'quantity.biomass': '지상부 바이오매스',
    'quantity.shannon': '수종 다양성',
    'quantity.conifer_fraction_sd': '침엽수 비율 — 불확실성',
    'quantity.biomass_sd': '지상부 바이오매스 — 불확실성',
    'quantity.shannon_sd': '수종 다양성 — 불확실성',
    'quantity.dominant_species_uncertainty': '우점 수종 — 불확실성',
    'quantity.plot_index': '할당된 표본점번호',
    'quantity.nearest_reference_distance': '최근접 참조 표본 거리',
    'place.south_korea': '대한민국',
    'place.jirisan': '지리산국립공원',
    'place.jeju': '제주도',
    'place.lower_geum': '금강 하류 유역',
    'mobile.settings_toggle': '지도 설정 표시 / 숨기기',
    'mobile.readout_toggle': '지점 상세정보 표시 / 숨기기',
    'mobile.info_toggle': '데이터셋 정보 표시 / 숨기기',
    'language.toggle': '한국어 / English',
    'basemap.control': '기본지도 설정',
    'basemap.map': '지도',
    'basemap.satellite': '위성',
    'basemap.labels': '라벨 표시',
    'location.button': '현재 위치 보기',
    'location.layer': '현재 위치',
    'legend.species_others': '기타 수종 (23종, 전체 면적의 5%)',
    'legend.plot_index': '색상은 할당된 표본점의 번호를 나타냅니다. 색상 표현에는 중복이 있습니다.',
    'range.pooled_p5_p95': '이 밴드를 포함하는 모든 레이어를 통합한 p5–p95',
    'range.definitional': '정의에 따른 범위',
    'unit.biomass': 'Mg ha⁻¹',
    'unit.biomass_sd': '앙상블 표준편차, Mg ha⁻¹',
    'unit.conifer_fraction': '흉고단면적 비율',
    'unit.conifer_fraction_sd': '앙상블 표준편차, 흉고단면적 비율',
    'unit.dominant_species_uncertainty': '1 − 투표 비율 차이',
    'unit.nearest_reference_distance': '1 − 말단 노드 근접도',
    'unit.shannon': 'Shannon H′',
    'unit.shannon_sd': '앙상블 표준편차, Shannon H′',
    'readout.prompt': '지도를 클릭하면 해당 셀의 모든 레이어 값을 확인할 수 있습니다.',
    'readout.loading': '불러오는 중…',
    'readout.no_data': '이 위치에는 데이터가 없습니다.',
    'readout.outside_extent': '산림 지도 범위 밖입니다.',
    'readout.unknown_species': '알 수 없는 수종',
    'readout.other_species': '주요 36종 외 기타 수종',
    'readout.conifer_fraction': '침엽수 비율',
    'readout.biomass': '지상부 바이오매스 (Mg ha⁻¹)',
    'readout.shannon': '수종 다양성 (Shannon H′)',
    'readout.dominant_species': '우점 수종',
    'readout.direct': '직접 예측',
    'readout.imputed': '표본점 할당',
    'readout.matched_plot': '할당된 참조 표본점',
    'readout.plot_id': '표본점 ID',
    'readout.distance': '거리',
    'coordinate.north': '북위',
    'coordinate.east': '동경',
    'note.spatial_pattern': '이 정보는 30m 해상도의 현장 실측값이 아닙니다. 광역적인 공간 ' +
      '패턴으로 해석하세요. 제작자는 현장 오류에 대한 모든 연락을 환영합니다.',
    'note.performance': '침엽수 비율과 우점 수종은 셀 단위로 읽을 수 있습니다. 바이오매스와 수종 ' +
      '다양성은 셀 단위 성능이 훨씬 낮으므로 1 km 이상으로 집계하여 읽도록 만들어졌습니다. 이 두 ' +
      '속성에서는 표본점 할당 방식의 성능이 매우 낮으므로 직접 예측 지도를 우선으로 사용하세요.',
    'note.uncertainty': '불확실성 레이어는 상대적인 신뢰도를 나타내며, 최근접 참조 거리는 해당 ' +
      '지점이 할당된 표본점과 얼마나 다른지를 나타냅니다. 둘 다 예측구간은 아닙니다.',
    'note.dominant_species': '우점 수종 정확도는 흔한 수종에 집중되어 있습니다. (희소한 수종일수록 ' +
      '정확도가 낮은 경향이 있습니다.) 수종별 정확도는 데이터셋의 ' +
      'dominant_species_table.csv에 포함되어 있습니다.',
    'note.plot_index': '표본점번호 할당 방식은 개별 속성의 정확도 일부를 포기하는 대신, 실제 ' +
      '조사 표본점에서 관측된 정보 간의 조합을 유지합니다. 각 지점에 할당된 표본점번호를 ' +
      '제8차 NFI 자료의 표본점 번호와 매칭하면 해당 조사구의 전체 조사 기록을 활용할 수 있습니다.'
  }
};

var currentLanguage = DEFAULT_LANGUAGE;
var localizedLabels = [];
var localizedButtons = [];

function tr(id) {
  var language = UI_TEXT[currentLanguage] || UI_TEXT.en;
  return language[id] === undefined ? UI_TEXT.en[id] : language[id];
}

function localizedLabel(id, style, url, prefix) {
  var item = {id: id, prefix: prefix || ''};
  item.widget = ui.Label(item.prefix + tr(id), style, url);
  localizedLabels.push(item);
  return item.widget;
}

function localizedButton(id, config) {
  config.label = tr(id);
  var item = {id: id, widget: ui.Button(config)};
  localizedButtons.push(item);
  return item.widget;
}

function refreshLocalizedWidgets() {
  localizedLabels.forEach(function (item) {
    item.widget.setValue(item.prefix + tr(item.id));
  });
  localizedButtons.forEach(function (item) {
    item.widget.setLabel(tr(item.id));
  });
}

function translatedOption(value, ids) {
  return ids[value] ? tr(ids[value]) : value;
}

function translatedOptions(values, ids) {
  return values.map(function (value) { return translatedOption(value, ids); });
}

function internalOption(display, values, ids) {
  for (var i = 0; i < values.length; i++) {
    if (translatedOption(values[i], ids) === display) { return values[i]; }
  }
  return display;
}

// The dataset's collection and plot-id table.
var COLLECTION = 'projects/ee-korea-forest-maps/assets/korea-forest-attributes-v2';
var PLOT_IDS = ee.FeatureCollection('projects/ee-korea-forest-maps/assets/plot_ids_v2');

var DOI = 'https://doi.org/10.5281/zenodo.23313068';
var CODE = 'https://github.com/kim-seunguk/korea-forest-attribute-maps';

// ---- what can be shown ------------------------------------------------------------------------
// A layer is a family and a quantity. The asset name follows from those two.
var FAMILY = {
  'Direct prediction': 'direct',
  'Plot-index imputation': 'imputation',
  'Basemap only': null
};

var QUANTITY = {
  direct: {
    'Coniferous fraction': 'conifer_fraction',
    'Dominant species': 'dominant_species',
    'Aboveground biomass': 'biomass',
    'Tree-species diversity': 'shannon',
    'Coniferous fraction — uncertainty': 'conifer_fraction_sd',
    'Aboveground biomass — uncertainty': 'biomass_sd',
    'Tree-species diversity — uncertainty': 'shannon_sd',
    'Dominant species — uncertainty': 'dominant_species_uncertainty'
  },
  imputation: {
    'Coniferous fraction': 'conifer_fraction',
    'Dominant species': 'dominant_species',
    'Aboveground biomass': 'biomass',
    'Tree-species diversity': 'shannon',
    'Matched plot index': 'plot_index',
    'Nearest-reference distance': 'nearest_reference_distance'
  }
};

var FAMILY_TEXT = {
  'Direct prediction': 'family.direct',
  'Plot-index imputation': 'family.imputation',
  'Basemap only': 'family.basemap'
};
var QUANTITY_TEXT = {
  'Coniferous fraction': 'quantity.conifer_fraction',
  'Dominant species': 'quantity.dominant_species',
  'Aboveground biomass': 'quantity.biomass',
  'Tree-species diversity': 'quantity.shannon',
  'Coniferous fraction — uncertainty': 'quantity.conifer_fraction_sd',
  'Aboveground biomass — uncertainty': 'quantity.biomass_sd',
  'Tree-species diversity — uncertainty': 'quantity.shannon_sd',
  'Dominant species — uncertainty': 'quantity.dominant_species_uncertainty',
  'Matched plot index': 'quantity.plot_index',
  'Nearest-reference distance': 'quantity.nearest_reference_distance'
};

var VIS_RANGE_TEXT = {
  'p5-p95 pooled over every layer carrying this band': 'range.pooled_p5_p95',
  'definitional': 'range.definitional'
};
var VIS_UNIT_TEXT = {
  'Mg ha-1': 'unit.biomass',
  'ensemble SD, Mg ha-1': 'unit.biomass_sd',
  'fraction of basal area': 'unit.conifer_fraction',
  'ensemble SD, fraction of basal area': 'unit.conifer_fraction_sd',
  '1 - vote margin': 'unit.dominant_species_uncertainty',
  '1 - terminal-node proximity': 'unit.nearest_reference_distance',
  'Shannon H′': 'unit.shannon',
  'ensemble SD, Shannon H′': 'unit.shannon_sd'
};

function translatedMetadata(value, ids) {
  return ids[value] ? tr(ids[value]) : value;
}

function visCaption(band) {
  var v = VIS[band];
  if (!v) { return ''; }
  return translatedMetadata(v.unit || '', VIS_UNIT_TEXT) + '   ' +
         translatedMetadata(v.range_from || '', VIS_RANGE_TEXT);
}

/** Asset id for a (family, band). */
function assetId(family, band) {
  var stem;
  if (band === 'plot_index') {
    stem = 'plot_index_30m';
  } else if (band === 'nearest_reference_distance') {
    stem = 'nearest_reference_distance_30m';
  } else if (family === 'imputation') {
    stem = 'imputed_' + band + '_30m';
  } else {
    stem = band + '_30m';
  }
  return COLLECTION + '/' + stem;
}

// Ramps read out of the ArcGIS project that drew the figures, so a value carries the same
// colour in the figures and here. Uncertainty and the nearest-reference distance share one:
// dark is less to be trusted.
var RAMPS = {
  conifer_fraction: ['1B7837', 'F7F7F7', '2166AC'],
  biomass: ['FEFE97', '13E22D', '00A2B8', '131A89'],
  shannon: ['000004', '1D1147', '51127C', '832682', 'B7377A', 'E65164', 'FC8761', 'FFC387', 'FDFEC0'],
  uncertainty: ['FCFBFD', 'EFEDF5', 'DADAEB', 'BCBDDC', '9E9AC8', '807DBA', '6A51A3', '54278F', '3F007D']
};
var MOSAIC = ['440154', '3B528B', '21908C', '5DC863', 'FDE725'];   // the plot index, whose value has no order

// Above 1 pulls colour into the low end of the range, below 1 into the high end. The figures use a
// histogram-equalised stretch, which has no equivalent here; this is the knob that comes closest.
var GAMMA = {};

function rampFor(band) {
  return isUncertainty(band) ? RAMPS.uncertainty : (RAMPS[band] || MOSAIC);
}

/** Layers where a high value is worse: they get the sequential ramp, not the viridis one. */
function isUncertainty(band) {
  return band.indexOf('_sd') > 0 || band.indexOf('uncertainty') > 0 ||
         band === 'nearest_reference_distance';
}

/** Rendered image plus the legend it needs. */
function styled(family, band) {
  var img = ee.Image(assetId(family, band));
  if (band === 'dominant_species') {
    return {image: img.visualize({min: 0, max: SPECIES.length - 1, palette: SPECIES}), kind: 'species'};
  }
  if (band === 'plot_index') {
    // The index carries no order, so it is coloured by its remainder mod 11, which shows where the
    // matched plot changes.
    return {image: img.mod(11).visualize({min: 0, max: 10, palette: MOSAIC.concat(MOSAIC)}),
            kind: 'mosaic'};
  }
  var v = VIS[band];
  var t = img.unitScale(v.min, v.max).clamp(0, 1);
  var g = GAMMA[band] || 1;
  if (g !== 1) { t = t.pow(1 / g); }
  return {image: t.visualize({min: 0, max: 1, palette: rampFor(band)}),
          kind: 'ramp', min: v.min, max: v.max, band: band};
}

// ---- maps -------------------------------------------------------------------------------------
// Both maps get the home view BEFORE they are linked. A Linker synchronises on change, so linking
// two maps that disagree lets whichever settles last win, and the one never given a centre carries
// the editor's default view of North America.
var HOME = {lon: 127.8, lat: 36.3, zoom: 7};
// The mainland, so it fills the map. Jeju sits far enough south that including it
// pushes the whole country into the middle third of the screen; it is one click away
// in the places list.
var HOME_BOX = ee.Geometry.Rectangle([126.0, 34.2, 129.6, 38.65], null, false);
var leftMap = ui.Map(), rightMap = ui.Map();
[leftMap, rightMap].forEach(function (m) {
  m.setOptions('SATELLITE');
  m.setControlVisibility({all: false, zoomControl: true, mapTypeControl: false});
  m.style().set({cursor: 'crosshair', stretch: 'both'});
  m.setCenter(HOME.lon, HOME.lat, HOME.zoom);
});
var mapLinker = ui.Map.Linker([leftMap, rightMap]);

// Keep the current view as client-side numbers. Re-parenting a map can reset its view, and using
// map.getCenter() during every layout change would require a server round trip.
var view = {lon: HOME.lon, lat: HOME.lat, zoom: HOME.zoom};
[leftMap, rightMap].forEach(function (m) {
  m.onChangeCenter(function (c) { view.lon = c.lon; view.lat = c.lat; });
  m.onChangeZoom(function (z) { view.zoom = z; });
});

function restoreView(map) {
  map.setCenter(view.lon, view.lat, view.zoom);
}

// ---- independent basemap controls -------------------------------------------------------------
var BASEMAP_NO_LABELS_ID = 'ROADMAP_NO_LABELS';
var BASEMAP_STYLES = {};
BASEMAP_STYLES[BASEMAP_NO_LABELS_ID] = [
  {elementType: 'labels', stylers: [{visibility: 'off'}]}
];
var BASEMAP_TYPES = ['map', 'satellite'];
var basemapControls = [];

function basemapTypeLabel(mapType) {
  return tr(mapType === 'map' ? 'basemap.map' : 'basemap.satellite');
}

function applyBasemap(control) {
  var mapTypeId = control.mapType === 'satellite' ?
    (control.labels ? 'HYBRID' : 'SATELLITE') :
    (control.labels ? 'ROADMAP' : BASEMAP_NO_LABELS_ID);
  control.map.setOptions({
    mapTypeId: mapTypeId,
    styles: BASEMAP_STYLES,
    types: ['ROADMAP', 'SATELLITE', 'HYBRID', BASEMAP_NO_LABELS_ID]
  });
}

function syncBasemapControlLanguage(control) {
  control.selector.items().reset(BASEMAP_TYPES.map(basemapTypeLabel));
  control.selector.setValue(basemapTypeLabel(control.mapType), false);
  control.labelsCheckbox.setLabel(tr('basemap.labels'));
}

function refreshBasemapControls() {
  basemapControls.forEach(syncBasemapControlLanguage);
}

function setBasemapMenuOpen(control, open) {
  basemapControls.forEach(function (item) {
    item.open = item === control ? open : false;
    item.body.style().set('shown', item.open);
  });
}

function makeBasemapControl(map, position, margin) {
  var control = {map: map, mapType: 'satellite', labels: false, open: false};
  control.selector = ui.Select({
    items: BASEMAP_TYPES.map(basemapTypeLabel),
    value: basemapTypeLabel(control.mapType),
    style: {stretch: 'horizontal', margin: '0 0 4px 0', fontSize: '12px'},
    onChange: function (display) {
      control.mapType = display === tr('basemap.map') ? 'map' : 'satellite';
      applyBasemap(control);
    }
  });
  control.labelsCheckbox = ui.Checkbox(
    tr('basemap.labels'), control.labels, function (checked) {
      control.labels = checked;
      applyBasemap(control);
    }, false, {stretch: 'horizontal', margin: '2px 0 0 0', fontSize: '12px'});
  control.body = ui.Panel([control.selector, control.labelsCheckbox],
                          ui.Panel.Layout.flow('vertical'),
                          {shown: false, width: '132px', padding: '6px', margin: '2px 0 0 0',
                           backgroundColor: '#FFFFFF', border: '1px solid #BDBDBD'});
  control.toggle = localizedButton('basemap.control', {
    imageUrl: BASEMAP_CONTROL_ICON,
    style: {width: '32px', height: '32px', padding: '0',
            margin: position === 'top-right' ? '0 0 0 auto' : '0'},
    onClick: function () {
      setBasemapMenuOpen(control, !control.open);
    }
  });
  control.panel = ui.Panel([control.toggle, control.body], ui.Panel.Layout.flow('vertical'),
                           {position: position, margin: margin, padding: '0'});
  basemapControls.push(control);
  map.add(control.panel);
  applyBasemap(control);
  return control;
}

// The wipe can cover the right edge of the left map, so its control sits beside the zoom buttons.
var leftBasemapControl = makeBasemapControl(leftMap, 'top-left', '8px 0 0 46px');
var rightBasemapControl = makeBasemapControl(rightMap, 'top-right', '8px');

// ---- current device location ------------------------------------------------------------------
var currentLocationPoint = null;
var locationButtons = [];

function appendCurrentLocationLayer(layers) {
  if (currentLocationPoint) {
    layers.push(ui.Map.Layer(currentLocationPoint, {color: 'E53935'}, tr('location.layer')));
  }
  return layers;
}

function refreshLocationButtons() {
  locationButtons.forEach(function (button) {
    button.setLabel(tr('location.button'));
  });
}

function findCurrentLocation(button) {
  button.setDisabled(true);
  ui.util.getCurrentPosition(function (point) {
    currentLocationPoint = point;
    sideStates.forEach(function (state) { state.refresh(); });
    var target = currentLayout === 'mobile' ? mobileSide.map :
      (desktopMapDisplay === 'Right map only' ? rightMap : leftMap);
    target.centerObject(point, CURRENT_LOCATION_ZOOM);
    button.setDisabled(false);
  }, function () {
    button.setDisabled(false);
  });
}

function makeLocationControl() {
  var button = ui.Button({
    label: tr('location.button'),
    imageUrl: CURRENT_LOCATION_ICON,
    style: {width: '28px', height: '28px', padding: '0', margin: '0'},
    onClick: findCurrentLocation
  });
  locationButtons.push(button);
  var panel = ui.Panel([button], ui.Panel.Layout.flow('horizontal'),
                       {position: 'bottom-right', margin: '6px', padding: '0'});
  return {button: button, panel: panel};
}

var leftLocationControl = makeLocationControl();
var rightLocationControl = makeLocationControl();
var leftLocationButton = leftLocationControl.button;
var rightLocationButton = rightLocationControl.button;
leftMap.add(leftLocationControl.panel);
rightMap.add(rightLocationControl.panel);

// ---- legend -----------------------------------------------------------------------------------
function rampLegend(spec, mobile) {
  var labelFont = mobile ? MOBILE_SMALL_FONT_SIZE : '12px';
  var bar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select('longitude'),
    params: {bbox: [0, 0, 1, 0.1], dimensions: '160x10', format: 'png',
             min: 0, max: 1, palette: rampFor(spec.band)},
    style: {stretch: 'horizontal', margin: '2px 0'}
  });
  return ui.Panel([
    bar,
    ui.Panel([ui.Label(String(spec.min), {fontSize: labelFont, margin: '0'}),
              ui.Label(String(spec.max), {fontSize: labelFont, margin: '0 0 0 auto'})],
             ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'})
  ], ui.Panel.Layout.flow('vertical'), mobile ? {stretch: 'horizontal'} : {width: '170px'});
}

function speciesLegend(mobile) {
  var labelFont = mobile ? MOBILE_BODY_FONT_SIZE : '13px';
  var rows = SPECIES_MAJOR.map(function (s) {
    var labelStyle = {fontSize: labelFont, margin: '1px 0'};
    if (currentLanguage === 'en') { labelStyle.fontStyle = 'italic'; }
    return ui.Panel([
      ui.Label('', {backgroundColor: '#' + s.hex, padding: '7px', margin: '1px 6px 1px 0'}),
      ui.Label(currentLanguage === 'ko' ? SPECIES_NAME_KO[s.code] : s.label, labelStyle)
    ], ui.Panel.Layout.flow('horizontal'));
  });
  rows.push(ui.Panel([
    ui.Label('', {backgroundColor: '#E1E1E1', padding: '7px', margin: '1px 6px 1px 0'}),
    ui.Label(tr('legend.species_others'), {fontSize: labelFont, margin: '1px 0'})
  ], ui.Panel.Layout.flow('horizontal')));
  return ui.Panel(rows, ui.Panel.Layout.flow('vertical'),
                  mobile ? {stretch: 'horizontal'} : {width: '244px'});
}

function legendWidget(spec, mobile) {
  return spec.kind === 'species' ? speciesLegend(mobile) :
    spec.kind === 'mosaic' ? ui.Label(tr('legend.plot_index'),
                                      {fontSize: mobile ? MOBILE_SMALL_FONT_SIZE : '12px',
                                       color: '#555'}) :
    rampLegend(spec, mobile);
}

// ---- one side ---------------------------------------------------------------------------------
function side(map, defaults) {
  // The state is independent of any one selector. Desktop has one selector set per map, while
  // mobile has one shared selector set that edits whichever map is active.
  var state = {
    map: map,
    family: defaults.family,
    quantity: defaults.quantity,
    wantedQuantity: defaults.quantity,
    spec: null,
    syncing: false
  };
  var familyNames = Object.keys(FAMILY);
  var familySel = ui.Select({
    items: translatedOptions(familyNames, FAMILY_TEXT),
    value: translatedOption(defaults.family, FAMILY_TEXT),
    style: {width: '240px'}
  });
  var quantSel = ui.Select({items: [], value: null, style: {width: '240px'}});
  var legendHolder = ui.Panel();
  var caption = ui.Label('', {fontSize: '12px', color: '#555', margin: '2px 0 0 0'});

  function quantityNames() {
    var family = FAMILY[state.family];
    return family ? Object.keys(QUANTITY[family]) : [];
  }

  function syncDesktopSelectors() {
    var names = quantityNames();
    state.syncing = true;
    familySel.items().reset(translatedOptions(familyNames, FAMILY_TEXT));
    familySel.setValue(translatedOption(state.family, FAMILY_TEXT), false);
    quantSel.items().reset(translatedOptions(names, QUANTITY_TEXT));
    if (names.length) {
      quantSel.setValue(translatedOption(state.quantity, QUANTITY_TEXT), false);
    }
    state.syncing = false;
  }

  function refresh() {
    var family = FAMILY[state.family];
    var band = family && QUANTITY[family][state.quantity];
    if (!band) {
      map.layers().reset(appendCurrentLocationLayer([]));
      legendHolder.widgets().reset([]);
      caption.setValue('');
      state.spec = null;
      return;
    }
    var spec = styled(family, band);
    state.spec = spec;
    map.layers().reset(appendCurrentLocationLayer([
      ui.Map.Layer(spec.image, {}, translatedOption(state.quantity, QUANTITY_TEXT))
    ]));
    legendHolder.widgets().reset([legendWidget(spec, false)]);
    caption.setValue(visCaption(band));
  }

  state.setFamily = function (label) {
    if (state.quantity) { state.wantedQuantity = state.quantity; }
    state.family = label;
    var names = quantityNames();
    if (names.length) {
      state.quantity = names.indexOf(state.wantedQuantity) >= 0 ? state.wantedQuantity : names[0];
    } else {
      state.quantity = null;
    }
    syncDesktopSelectors();
    refresh();
  };

  state.setQuantity = function (label) {
    if (!label) { return; }
    state.quantity = label;
    state.wantedQuantity = label;
    syncDesktopSelectors();
    refresh();
  };

  state.quantityNames = quantityNames;
  state.refresh = refresh;
  state.updateLanguage = function () {
    syncDesktopSelectors();
    refresh();
  };

  familySel.onChange(function (label) {
    if (!state.syncing) { state.setFamily(internalOption(label, familyNames, FAMILY_TEXT)); }
  });
  quantSel.onChange(function (label) {
    if (!state.syncing) {
      state.setQuantity(internalOption(label, quantityNames(), QUANTITY_TEXT));
    }
  });

  var title = localizedLabel(defaults.titleKey,
                             {fontWeight: 'bold', fontSize: '13px', margin: '0 0 3px 0'});
  state.panel = ui.Panel([title, familySel, quantSel, legendHolder, caption],
                         ui.Panel.Layout.flow('vertical'),
                         {margin: '0 0 16px 0', width: '278px'});
  syncDesktopSelectors();
  refresh();
  return state;
}

// ---- click: read every layer at one cell, in one request --------------------------------------
/** Every published band, named family_band, as a single image. */
function fullStack() {
  var bands = [];
  ['direct', 'imputation'].forEach(function (family) {
    Object.keys(QUANTITY[family]).forEach(function (label) {
      var band = QUANTITY[family][label];
      bands.push(ee.Image(assetId(family, band)).rename(family + '_' + band));
    });
  });
  return ee.Image.cat(bands);
}
var STACK = fullStack();

function makeReadout(mobile) {
  return ui.Panel([ui.Label(tr('readout.prompt'),
                            {fontSize: mobile ? MOBILE_BODY_FONT_SIZE : '13px', color: '#555'})],
                  ui.Panel.Layout.flow('vertical'),
                  mobile ? {margin: '0', stretch: 'horizontal'} :
                           {margin: '0 0 16px 0', width: '278px'});
}

var desktopReadout = makeReadout(false);
var mobileReadout = makeReadout(true);
var READOUTS = [desktopReadout, mobileReadout];
var lastReadoutPoint = null;

function resetReadouts(makeWidgets) {
  READOUTS.forEach(function (panel) {
    panel.widgets().reset(makeWidgets(panel === mobileReadout));
  });
}

function showReadoutPrompt() {
  resetReadouts(function (mobile) {
    return [ui.Label(tr('readout.prompt'),
                     {fontSize: mobile ? MOBILE_BODY_FONT_SIZE : '13px', color: '#555'})];
  });
}

function fmt(v, digits) {
  return (v === null || v === undefined) ? '—' : Number(v).toFixed(digits);
}

/** One readout value in the form readoutCard renders: the text, italics for a scientific name, and
 *  a highlight on the family whose value is the one to read for that attribute. */
function readoutValue(part, highlight) {
  return {text: part.text,
          parts: [{text: part.text, italic: part.italic || false, highlight: highlight || false}]};
}

function numericValue(value, show, highlight) {
  return readoutValue({text: show(value), italic: false}, highlight);
}

function readoutCard(title, items, mobile) {
  var bodyFont = mobile ? MOBILE_BODY_FONT_SIZE : '13px';
  var keyWidth = mobile ? '88px' : '72px';
  var widgets = [ui.Label(title, {
    fontSize: bodyFont, fontWeight: 'bold', margin: '0 0 4px 0'
  })];
  var labels = {};
  items.forEach(function (item) {
    var valueLabel;
    if (item.parts) {
      var partLabels = item.parts.map(function (part) {
        var style = {
          fontSize: bodyFont, margin: part.separator ? '0 4px' : '0'
        };
        if (part.italic) { style.fontStyle = 'italic'; }
        if (part.highlight) {
          style.fontWeight = 'bold';
          style.color = READOUT_HIGHLIGHT_COLOR;
        }
        return ui.Label(part.text, style);
      });
      valueLabel = ui.Panel(partLabels, ui.Panel.Layout.flow('horizontal', true),
                            {stretch: 'horizontal'});
    } else {
      valueLabel = ui.Label(item.value, {
        fontSize: bodyFont, margin: '0', stretch: 'horizontal'
      });
    }
    if (item.id) { labels[item.id] = valueLabel; }
    widgets.push(ui.Panel([
      ui.Label(item.key, {fontSize: bodyFont, color: '#555', width: keyWidth, margin: '0'}),
      valueLabel
    ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));
  });
  return {
    panel: ui.Panel(widgets, ui.Panel.Layout.flow('vertical'), {
      stretch: 'horizontal', padding: mobile ? '8px' : '6px', margin: '0 0 6px 0',
      backgroundColor: '#F5F5F5', border: '1px solid #DDDDDD'
    }),
    labels: labels
  };
}

function speciesReadoutPart(value) {
  if (value === null || value === undefined) { return {text: '—', italic: false}; }
  var number = Number(value);
  var code = Math.round(number);
  if (!isFinite(number) || Math.abs(number - code) > 1e-6 || code < 0 || code > 36) {
    return {text: tr('readout.unknown_species'), italic: false};
  }
  if (code === 36) { return {text: tr('readout.other_species'), italic: false}; }
  return currentLanguage === 'ko' ?
    {text: SPECIES_NAME_KO[code], italic: false} :
    {text: SPECIES_NAME[code], italic: true};
}

function speciesValue(value, highlight) {
  return readoutValue(speciesReadoutPart(value), highlight);
}

function onClick(pt) {
  lastReadoutPoint = {lat: pt.lat, lon: pt.lon};
  if (currentLayout === 'mobile') { setMobileReadoutOpen(true); }
  resetReadouts(function (mobile) {
    return [ui.Label(tr('readout.loading'),
                     {fontSize: mobile ? MOBILE_BODY_FONT_SIZE : '13px'})];
  });
  var geom = ee.Geometry.Point([pt.lon, pt.lat]);
  STACK.reduceRegion({reducer: ee.Reducer.first(), geometry: geom, scale: 30})
    .evaluate(function (v) {
      if (!v) {
        resetReadouts(function (mobile) {
          return [ui.Label(tr('readout.no_data'),
                           {fontSize: mobile ? MOBILE_BODY_FONT_SIZE : '13px'})];
        });
        return;
      }

      var any = false;
      var attributeValues = [];
      var attributeLabels = {
        conifer_fraction: 'readout.conifer_fraction',
        biomass: 'readout.biomass',
        shannon: 'readout.shannon',
        dominant_species: 'readout.dominant_species'
      };
      // One value per attribute is highlighted: the direct prediction for the three continuous
      // attributes, the matched plot's species for dominant species.
      ['conifer_fraction', 'biomass', 'shannon', 'dominant_species'].forEach(function (band) {
        var species = band === 'dominant_species';
        var d, i;
        if (species) {
          d = speciesValue(v['direct_' + band], false);
          i = speciesValue(v['imputation_' + band], true);
        } else {
          var show = function (x) { return fmt(x, band === 'biomass' ? 1 : 3); };
          d = numericValue(v['direct_' + band], show, true);
          i = numericValue(v['imputation_' + band], show, false);
        }
        if (d.text === '\u2014' && i.text === '\u2014') { return; }
        any = true;
        attributeValues.push({title: tr(attributeLabels[band]), direct: d, imputed: i});
      });

      var dist = fmt(v['imputation_nearest_reference_distance'], 3);
      var plotIndex = v['imputation_plot_index'];
      var hasPlot = plotIndex !== null && plotIndex !== undefined;
      var plotLabels = [];

      READOUTS.forEach(function (panel) {
        var mobile = panel === mobileReadout;
        var bodyFont = mobile ? MOBILE_BODY_FONT_SIZE : '13px';
        var coordinates = currentLanguage === 'ko' ?
          tr('coordinate.north') + ' ' + fmt(pt.lat, 5) + '°,  ' +
          tr('coordinate.east') + ' ' + fmt(pt.lon, 5) + '°' :
          fmt(pt.lat, 5) + ' N,  ' + fmt(pt.lon, 5) + ' E';
        var rows = [ui.Label(coordinates,
                             {fontSize: mobile ? MOBILE_TITLE_FONT_SIZE : '14px',
                              fontWeight: 'bold', margin: '0 0 8px 0'})];
        attributeValues.forEach(function (attribute) {
          rows.push(readoutCard(attribute.title, [
            {key: tr('readout.direct'), value: attribute.direct.text, parts: attribute.direct.parts},
            {key: tr('readout.imputed'), value: attribute.imputed.text, parts: attribute.imputed.parts}
          ], mobile).panel);
        });
        if (hasPlot) {
          var plotCard = readoutCard(tr('readout.matched_plot'), [
            {key: tr('readout.plot_id'), value: tr('readout.loading'), id: 'plotId'},
            {key: tr('readout.distance'), value: dist}
          ], mobile);
          rows.push(plotCard.panel);
          plotLabels.push(plotCard.labels.plotId);
        }
        if (!any) {
          rows.push(ui.Label(tr('readout.outside_extent'), {fontSize: bodyFont}));
        }
        panel.widgets().reset(rows);
      });

      if (hasPlot) {
        var id = ee.String(PLOT_IDS.filter(ee.Filter.eq('plot_index', Math.round(plotIndex)))
                           .first().get('nfi_plot_id'));
        id.evaluate(function (got) {
          var name = (got === null || got === undefined) ? '\u2014' : String(got);
          plotLabels.forEach(function (lab) {
            lab.setValue(name);
          });
        });
      }
    });
}

leftMap.onClick(onClick);
rightMap.onClick(onClick);

// ---- places ------------------------------------------------------------------------------------
// Bounds read from the zoom rectangles the map figures use, in EPSG:5179.
var PLACES = [
  {name: 'South Korea', bounds: null},
  {name: 'Jirisan National Park', bounds: [995689.2, 1690744.1, 1029510.7, 1716615.9]},
  {name: 'Jeju Island', bounds: [869275.0, 1457040.1, 950275.0, 1518999.9]},
  {name: 'Lower Geum River basin', bounds: [926491.0, 1774630.6, 939016.0, 1784211.4]}
];
var PLACE_TEXT = {
  'South Korea': 'place.south_korea',
  'Jirisan National Park': 'place.jirisan',
  'Jeju Island': 'place.jeju',
  'Lower Geum River basin': 'place.lower_geum'
};

var placeState = 'South Korea';
var placeSelectors = [];

function placeNames() {
  return PLACES.map(function (p) { return p.name; });
}

function syncPlaceSelectors() {
  var names = placeNames();
  placeSelectors.forEach(function (selector) {
    selector.items().reset(translatedOptions(names, PLACE_TEXT));
    selector.setValue(translatedOption(placeState, PLACE_TEXT), false);
  });
}

function goToPlace(name, source) {
  placeState = name;
  placeSelectors.forEach(function (selector) {
    if (selector !== source) {
      selector.setValue(translatedOption(name, PLACE_TEXT), false);
    }
  });
  var p = PLACES.filter(function (x) { return x.name === name; })[0];
  var target = currentLayout === 'mobile' && mobileSide ? mobileSide.map : leftMap;
  if (!p.bounds) { target.centerObject(HOME_BOX); return; }
  var box = ee.Geometry.Rectangle(p.bounds, 'EPSG:5179', false);
  target.centerObject(box);
}

function makePlaceSelector(mobile) {
  var names = placeNames();
  var selector = ui.Select({
    items: translatedOptions(names, PLACE_TEXT),
    value: translatedOption(placeState, PLACE_TEXT),
    style: mobile ? {stretch: 'horizontal', fontSize: MOBILE_BODY_FONT_SIZE} :
                    {stretch: 'horizontal'}
  });
  selector.onChange(function (name) {
    goToPlace(internalOption(name, names, PLACE_TEXT), selector);
  });
  placeSelectors.push(selector);
  return selector;
}

// ---- one map or two ----------------------------------------------------------------------------
var splitPanel = ui.SplitPanel({
  firstPanel: ui.Panel(), secondPanel: ui.Panel(),
  orientation: 'horizontal', wipe: true, style: {stretch: 'both'}
});
var mapHolder = ui.Panel([], ui.Panel.Layout.flow('vertical'), {stretch: 'both'});
var desktopMapDisplay = 'Compare two maps';
var DESKTOP_LAYOUTS = ['Compare two maps', 'Left map only', 'Right map only'];
var LAYOUT_TEXT = {
  'Compare two maps': 'layout.compare',
  'Left map only': 'layout.left_only',
  'Right map only': 'layout.right_only'
};

function releaseMapsFromDesktop() {
  mapHolder.widgets().reset([]);
  splitPanel.setFirstPanel(ui.Panel());
  splitPanel.setSecondPanel(ui.Panel());
}

function showDesktopMaps() {
  var saved = {lon: view.lon, lat: view.lat, zoom: view.zoom};
  releaseMapsFromDesktop();
  var shown = desktopMapDisplay === 'Right map only' ? rightMap : leftMap;
  if (desktopMapDisplay === 'Compare two maps') {
    splitPanel.setFirstPanel(leftMap);
    splitPanel.setSecondPanel(rightMap);
    mapHolder.widgets().reset([splitPanel]);
  } else {
    mapHolder.widgets().reset([shown]);
  }
  view.lon = saved.lon;
  view.lat = saved.lat;
  view.zoom = saved.zoom;
  restoreView(shown);
}

var layoutSel = ui.Select({
  items: translatedOptions(DESKTOP_LAYOUTS, LAYOUT_TEXT),
  value: translatedOption(desktopMapDisplay, LAYOUT_TEXT),
  style: {stretch: 'horizontal'},
  onChange: function (v) {
    desktopMapDisplay = internalOption(v, DESKTOP_LAYOUTS, LAYOUT_TEXT);
    if (currentLayout === 'desktop') { showDesktopMaps(); }
  }
});

function syncLayoutSelector() {
  layoutSel.items().reset(translatedOptions(DESKTOP_LAYOUTS, LAYOUT_TEXT));
  layoutSel.setValue(translatedOption(desktopMapDisplay, LAYOUT_TEXT), false);
}

// ---- usage notes -------------------------------------------------------------------------------
var NOTE_IDS = [
  'note.spatial_pattern',
  'note.performance',
  'note.uncertainty',
  'note.dominant_species',
  'note.plot_index'
];

function desktopButtonRow(button, margin) {
  return ui.Panel([button], ui.Panel.Layout.flow('horizontal'),
                  {stretch: 'horizontal', margin: margin});
}

function makeLanguageButton(mobile) {
  var button = ui.Button({
    label: tr('language.toggle'),
    style: mobile ? {stretch: 'horizontal', fontSize: MOBILE_SMALL_FONT_SIZE, padding: '6px',
                     margin: '6px 0 0 0'} :
                    {width: '200px', margin: '0'},
    onClick: toggleLanguage
  });
  return mobile ? button : desktopButtonRow(button, '4px 0 6px 0');
}

function makeInfoPanel(mobile) {
  var titleFont = mobile ? MOBILE_TITLE_FONT_SIZE : '16px';
  var bodyFont = mobile ? MOBILE_BODY_FONT_SIZE : '13px';
  var smallFont = mobile ? MOBILE_SMALL_FONT_SIZE : '12px';
  var description = localizedLabel(
    'info.description',
    {fontSize: bodyFont, margin: '2px 0 2px 0'}
  );
  var widgets = mobile ? [description] : [
    localizedLabel('header.title', {fontWeight: 'bold', fontSize: titleFont, margin: '0'}),
    description,
    localizedLabel('info.author',
                   {fontSize: smallFont, color: '#555', margin: '0 0 6px 0'})
  ];
  if (!mobile) {
    widgets.push(makeLanguageButton(false));
    widgets.push(layoutSel);
  }
  widgets.push(makePlaceSelector(mobile));
  widgets.push(localizedLabel('info.data_link',
                              {fontSize: bodyFont, margin: '6px 0 0 0'}, DOI));
  widgets.push(localizedLabel('info.code_link', {fontSize: bodyFont, margin: '0'}, CODE));

  var panel = ui.Panel(widgets, ui.Panel.Layout.flow('vertical'),
                       mobile ? {stretch: 'horizontal'} :
                                {margin: '0 0 16px 0', width: '278px'});
  var notes = ui.Panel([], ui.Panel.Layout.flow('vertical'), {shown: false, margin: '0'});
  NOTE_IDS.forEach(function (id) {
    notes.add(localizedLabel(id,
                             {fontSize: smallFont, color: '#444', margin: '3px 0 0 0'},
                             null, '• '));
  });
  var notesButton = localizedButton('info.notes_toggle', {
    style: mobile ? {margin: '6px 0 0 0', stretch: 'horizontal',
                     fontSize: MOBILE_BODY_FONT_SIZE, padding: '8px'} :
                    {width: '200px', margin: '0'},
    onClick: function () { notes.style().set('shown', !notes.style().get('shown')); }
  });
  panel.add(mobile ? notesButton : desktopButtonRow(notesButton, '6px 0 0 0'));
  panel.add(notes);
  return panel;
}

var desktopInfo = makeInfoPanel(false);
var mobileInfo = makeInfoPanel(true);


// ---- responsive layout -------------------------------------------------------------------------
// Desktop keeps the original fixed sidebars and wipe. Mobile uses separate controls and content
// panels, so no selector, legend or readout widget ever has two parents.
var leftSide = side(leftMap, {family: 'Direct prediction', quantity: 'Aboveground biomass',
                              titleKey: 'map.left_title'});
var rightSide = side(rightMap, {family: 'Plot-index imputation', quantity: 'Dominant species',
                                titleKey: 'map.right_title'});
var sideStates = [leftSide, rightSide];
var mobileSide = rightSide;

function bar(widgets) {
  return ui.Panel(widgets, ui.Panel.Layout.flow('vertical'),
                  {width: '300px', padding: '8px', backgroundColor: '#FAFAFA'});
}

var desktopRoot = ui.Panel([
  bar([leftSide.panel, rightSide.panel]),
  mapHolder,
  bar([desktopInfo, desktopReadout])
], ui.Panel.Layout.flow('horizontal'), {stretch: 'both'});

// Mobile exposes one map only. Its controls edit the desktop right-map state; the
// desktop left-map state remains intact and returns when the desktop wipe is shown again.
var mobileControlSync = false;
var mobileMapHolder = ui.Panel([], ui.Panel.Layout.flow('vertical'),
                               {stretch: 'horizontal', width: '100%', maxWidth: '100%',
                                minWidth: '0px', height: '420px'});
var mobileLegendHolder = ui.Panel([], ui.Panel.Layout.flow('vertical'), {stretch: 'horizontal'});
var mobileCaption = ui.Label('',
                             {fontSize: MOBILE_SMALL_FONT_SIZE, color: '#555', margin: '4px 0 0 0'});

var mobileFamilySel = ui.Select({
  items: translatedOptions(Object.keys(FAMILY), FAMILY_TEXT),
  style: {stretch: 'horizontal', fontSize: MOBILE_BODY_FONT_SIZE}
});
var mobileQuantSel = ui.Select({
  items: [], style: {stretch: 'horizontal', fontSize: MOBILE_BODY_FONT_SIZE}
});

function renderMobileLegend() {
  var state = mobileSide;
  mobileLegendHolder.widgets().reset(state.spec ? [legendWidget(state.spec, true)] : []);
  var family = FAMILY[state.family];
  var band = family && state.quantity ? QUANTITY[family][state.quantity] : null;
  mobileCaption.setValue(band ? visCaption(band) : '');
}

function syncMobileControls() {
  var state = mobileSide;
  var names = state.quantityNames();
  mobileControlSync = true;
  mobileFamilySel.items().reset(translatedOptions(Object.keys(FAMILY), FAMILY_TEXT));
  mobileFamilySel.setValue(translatedOption(state.family, FAMILY_TEXT), false);
  mobileQuantSel.items().reset(translatedOptions(names, QUANTITY_TEXT));
  if (names.length) {
    mobileQuantSel.setValue(translatedOption(state.quantity, QUANTITY_TEXT), false);
  }
  mobileControlSync = false;
  renderMobileLegend();
}

mobileFamilySel.onChange(function (label) {
  if (!mobileControlSync) {
    mobileSide.setFamily(internalOption(label, Object.keys(FAMILY), FAMILY_TEXT));
    syncMobileControls();
  }
});
mobileQuantSel.onChange(function (label) {
  if (!mobileControlSync) {
    mobileSide.setQuantity(internalOption(label, mobileSide.quantityNames(), QUANTITY_TEXT));
    syncMobileControls();
  }
});

var mobileSettingsBody = ui.Panel([
  localizedLabel('settings.family',
                 {fontSize: MOBILE_SMALL_FONT_SIZE, color: '#555', margin: '0'}),
  mobileFamilySel,
  localizedLabel('settings.quantity',
                 {fontSize: MOBILE_SMALL_FONT_SIZE, color: '#555', margin: '8px 0 0 0'}),
  mobileQuantSel,
  mobileCaption
], ui.Panel.Layout.flow('vertical'),
   {shown: false, stretch: 'horizontal', padding: '8px', backgroundColor: '#FAFAFA'});

var mobileLegendBody = ui.Panel([mobileLegendHolder], ui.Panel.Layout.flow('vertical'),
                                {shown: true, stretch: 'horizontal', margin: '10px 0 0 0'});
mobileSettingsBody.add(mobileLegendBody);

var mobileSettingsToggle = localizedButton('mobile.settings_toggle', {
  style: {stretch: 'horizontal', margin: '0 0 6px 0', fontSize: MOBILE_BODY_FONT_SIZE,
          padding: '8px'},
  onClick: function () {
    mobileSettingsBody.style().set('shown', !mobileSettingsBody.style().get('shown'));
  }
});

var mobileReadoutBody = ui.Panel([mobileReadout], ui.Panel.Layout.flow('vertical'),
                                 {shown: false, stretch: 'horizontal', padding: '8px'});
function setMobileReadoutOpen(open) {
  mobileReadoutBody.style().set('shown', open);
}
var mobileReadoutToggle = localizedButton('mobile.readout_toggle', {
  style: {stretch: 'horizontal', margin: '0 0 6px 0', fontSize: MOBILE_BODY_FONT_SIZE,
          padding: '8px'},
  onClick: function () {
    setMobileReadoutOpen(!mobileReadoutBody.style().get('shown'));
  }
});

var mobileInfoBody = ui.Panel([mobileInfo], ui.Panel.Layout.flow('vertical'),
                              {shown: false, stretch: 'horizontal', padding: '8px'});
var mobileInfoToggle = localizedButton('mobile.info_toggle', {
  style: {stretch: 'horizontal', margin: '0', fontSize: MOBILE_BODY_FONT_SIZE, padding: '8px'},
  onClick: function () {
    mobileInfoBody.style().set('shown', !mobileInfoBody.style().get('shown'));
  }
});

var mobileLanguageButton = makeLanguageButton(true);
var mobileHeader = ui.Panel([
  localizedLabel('header.title',
                 {fontWeight: 'bold', fontSize: MOBILE_TITLE_FONT_SIZE, margin: '0'}),
  localizedLabel('header.author_short',
                 {fontSize: MOBILE_SMALL_FONT_SIZE, color: '#555', margin: '2px 0 0 0'}),
  mobileLanguageButton
], ui.Panel.Layout.flow('vertical'),
   {stretch: 'horizontal', padding: '10px', backgroundColor: '#FAFAFA'});

var mobileActionPanel = ui.Panel([
  mobileSettingsToggle,
  mobileReadoutToggle,
  mobileInfoToggle
], ui.Panel.Layout.flow('vertical'),
   {stretch: 'horizontal', padding: '0 10px 10px 10px', backgroundColor: '#FAFAFA'});

var mobileRoot = ui.Panel([
  mobileHeader,
  mobileActionPanel,
  mobileSettingsBody,
  mobileReadoutBody,
  mobileInfoBody,
  mobileMapHolder
], ui.Panel.Layout.flow('vertical'),
   {stretch: 'both', width: '100%', maxWidth: '100%', minWidth: '0px'});

function setLanguage(language) {
  if (language !== 'en' && language !== 'ko') { return; }
  if (language === currentLanguage) { return; }
  currentLanguage = language;
  refreshLocalizedWidgets();
  refreshBasemapControls();
  refreshLocationButtons();
  syncLayoutSelector();
  syncPlaceSelectors();
  sideStates.forEach(function (state) { state.updateLanguage(); });
  syncMobileControls();
  if (lastReadoutPoint) {
    onClick(lastReadoutPoint);
  } else {
    showReadoutPrompt();
  }
}

function toggleLanguage() {
  setLanguage(currentLanguage === 'en' ? 'ko' : 'en');
}

function useMobileLayout(device) {
  if (LAYOUT_TEST_MODE === 'mobile') { return true; }
  if (LAYOUT_TEST_MODE === 'desktop') { return false; }
  return device.is_mobile ||
         device.width <= MOBILE_BREAKPOINT ||
         (device.is_tablet && device.is_portrait);
}

function sizeMapForResponsiveContainer(map) {
  // Explicit percentage sizing plus a zero minimum lets ui.Map shrink after it has previously
  // occupied a wide desktop SplitPanel. Stretch alone can retain that earlier intrinsic width.
  map.style().set({stretch: 'both', width: '100%', maxWidth: '100%', minWidth: '0px',
                   height: '100%'});
}

function updateMobileMapSize(device) {
  // Code Editor and device emulation can report the two viewport axes differently after rotation.
  // In portrait use the long edge; in landscape use the short edge. This also makes the map occupy
  // most of the first portrait screen instead of leaving a large unused area below it.
  var portrait = device.is_portrait && !device.is_landscape;
  var viewportHeight = portrait ? Math.max(device.width, device.height) :
                                  Math.min(device.width, device.height);
  var ratio = portrait ? MOBILE_MAP_PORTRAIT_RATIO : MOBILE_MAP_LANDSCAPE_RATIO;
  var minimum = portrait ? 420 : 280;
  var height = Math.max(minimum, Math.floor(viewportHeight * ratio));
  mobileMapHolder.style().set({height: height + 'px', width: '100%', maxWidth: '100%',
                               minWidth: '0px', stretch: 'horizontal'});
  [leftMap, rightMap].forEach(sizeMapForResponsiveContainer);
}

var currentLayout = null;
function applyResponsiveLayout(device) {
  var next = useMobileLayout(device) ? 'mobile' : 'desktop';
  if (next === currentLayout) {
    if (next === 'mobile') { updateMobileMapSize(device); }
    return;
  }

  var saved = {lon: view.lon, lat: view.lat, zoom: view.zoom};
  ui.root.widgets().reset([]);
  mobileMapHolder.widgets().reset([]);
  releaseMapsFromDesktop();
  view.lon = saved.lon;
  view.lat = saved.lat;
  view.zoom = saved.zoom;
  currentLayout = next;

  if (next === 'mobile') {
    updateMobileMapSize(device);
    mobileMapHolder.widgets().reset([mobileSide.map]);
    syncMobileControls();
    ui.root.widgets().reset([mobileRoot]);
    restoreView(mobileSide.map);
  } else {
    showDesktopMaps();
    ui.root.widgets().reset([desktopRoot]);
    restoreView(desktopMapDisplay === 'Right map only' ? rightMap : leftMap);
  }
}

// Initialise the original home extent, then let onResize choose the first layout and every later
// layout transition. The callback also runs once when the script starts.
leftMap.centerObject(HOME_BOX);
ui.root.onResize(applyResponsiveLayout);
