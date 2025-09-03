import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/Geomorphometry.jl/dev/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Home', link: '/index' },
{ text: 'Getting started', collapsed: false, items: [
{ text: 'Installation', link: '/installation' },
{ text: 'Usage', link: '/usage' },
{ text: 'Experimental', link: '/experimental' }]
 },
{ text: 'Background', collapsed: false, items: [
{ text: 'Concepts', link: '/concepts' },
{ text: 'Future plans', link: '/todo' }]
 },
{ text: 'Reference', collapsed: false, items: [
{ text: 'Validation', link: '/validation' },
{ text: 'API', link: '/reference' },
{ text: 'Changelog', link: '/CHANGELOG' },
{ text: 'Bibliography', link: '/bibliography' }]
 }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/Geomorphometry.jl/dev/',// TODO: replace this in makedocs!
  title: 'Geomorphometry.jl',
  description: 'Documentation for Geomorphometry.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'favicon.ico' }],
    ['script', { src: `${getBaseRepository(baseTemp.base)}versions.js` }],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', { src: `${baseTemp.base}siteinfo.js` }]
  ],
  ignoreDeadLinks: true,

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
        md.use(mathjax3),
        md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    }
  },
  themeConfig: {
    outline: 'deep',
    logo: '/logo.svg',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Getting started', collapsed: false, items: [
{ text: 'Installation', link: '/installation' },
{ text: 'Usage', link: '/usage' },
{ text: 'Experimental', link: '/experimental' }]
 },
{ text: 'Background', collapsed: false, items: [
{ text: 'Concepts', link: '/concepts' },
{ text: 'Future plans', link: '/todo' }]
 },
{ text: 'Reference', collapsed: false, items: [
{ text: 'Validation', link: '/validation' },
{ text: 'API', link: '/reference' },
{ text: 'Changelog', link: '/CHANGELOG' },
{ text: 'Bibliography', link: '/bibliography' }]
 }
]
,
    editLink: { pattern: "https://github.com/Deltares/Geomorphometry.jl/edit/main/docs/src/:path" },
    socialLinks: [
      // { icon: 'github', link: 'https://github.com/Deltares/Geomorphometry.jl' },
      { icon: 'linkedin', link: 'https://www.linkedin.com/in/mjpronk/' },
      { icon: 'mastodon', link: 'https://fosstodon.org/@evetion' },
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
