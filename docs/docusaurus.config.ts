import { themes as prismThemes } from 'prism-react-renderer';
import type { Config } from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';

const config: Config = {
  title: 'opentfraw-reader reference',
  tagline: 'Rust and Python reader for Thermo Fisher RAW mass spectrometry files',

  markdown: {
    mermaid: true,
    hooks: {
      onBrokenMarkdownLinks: 'warn',
    },
  },
  plugins: ['docusaurus-plugin-llms-txt'],
  themes: ['@docusaurus/theme-mermaid'],

  url: 'http://localhost:25812',
  baseUrl: '/',

  organizationName: 'pgarrett-scripps',
  projectName: 'opentfraw',

  onBrokenLinks: 'throw',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      {
        docs: {
          routeBasePath: '/',
          sidebarPath: './sidebars.ts',
          editUrl: 'https://github.com/pgarrett-scripps/opentfraw/tree/main/docs/',
        },
        blog: false,
        sitemap: {
          changefreq: 'weekly',
          priority: 0.5,
          filename: 'sitemap.xml',
        },
        theme: {
          customCss: './src/css/custom.css',
        },
      } satisfies Preset.Options,
    ],
  ],

  themeConfig: {
    metadata: [
      { name: 'keywords', content: 'OpenTFRaw, Thermo Fisher, RAW, mass spectrometry, Orbitrap, Rust, Python' },
      { name: 'description', content: 'OpenTFRaw is a Rust and Python reader for Thermo Fisher RAW mass spectrometry files.' },
    ],
    colorMode: {
      defaultMode: 'dark',
      disableSwitch: false,
      respectPrefersColorScheme: true,
    },
    announcementBar: {
      id: 'archived-reference',
      content: 'Inherited format reference. See the repository README and REVIEW.md for qualified fork behavior. This documentation is served locally only.',
      isCloseable: false,
    },
    navbar: {
      title: 'opentfraw-reader reference',
      items: [
        { href: 'https://github.com/pgarrett-scripps/opentfraw', label: 'Fork source', position: 'right' },
        { href: 'https://github.com/Sigilweaver/OpenTFRaw', label: 'Upstream source', position: 'right' },
      ],
    },
    footer: {
      style: 'dark',
      copyright: 'Original OpenTFRaw material copyright 2026 Sigilweaver Holdings LLC. Fork changes by opentfraw-reader contributors. Code Apache-2.0. Inherited documentation retains CC-BY-SA 4.0 attribution.',
    },
    prism: {
      theme: prismThemes.github,
      darkTheme: prismThemes.dracula,
      additionalLanguages: ['rust', 'toml', 'bash'],
    },
  } satisfies Preset.ThemeConfig,
};

export default config;
